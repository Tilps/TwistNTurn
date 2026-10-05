using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using TwistNTurn;

namespace TwistNTurn.Web.Services
{
    public class GameService
    {
        public Mesh? CurrentMesh { get; private set; }
        public UndoTree UndoTree { get; private set; } = new();
        public HashSet<int> MarkedEdges { get; private set; } = new();

        public MeshType CurrentType { get; set; } = MeshType.Square;
        public string SizeText { get; set; } = "10x10";
        public int Difficulty { get; set; } = 1; // Easy

        public bool DisallowFalseMove { get; set; } = false;
        public bool ShowCellColors { get; set; } = false;

        public bool IsGenerating { get; private set; }
        public int PrunedProgress { get; private set; }
        public int TotalIntersections { get; private set; }

        public DateTime? StartTime { get; private set; }
        public TimeSpan ElapsedTime => StartTime.HasValue ? DateTime.UtcNow - StartTime.Value : TimeSpan.Zero;
        public bool IsSolved { get; private set; }

        public event Action? OnStateChanged;

        public bool CanUndo => UndoTree.CanUndo;
        public bool CanRedo => UndoTree.CanRedo;

        public void NotifyChanged()
        {
            OnStateChanged?.Invoke();
        }

        private Mesh? generatingMesh;
        private CancellationTokenSource? generationCts;
        public DateTime? GenerationStartTime { get; private set; }

        public double GenerationElapsedSeconds => GenerationStartTime.HasValue
            ? (DateTime.UtcNow - GenerationStartTime.Value).TotalSeconds
            : 0;

        public void CancelGeneration()
        {
            if (generatingMesh != null)
            {
                generatingMesh.AbortPrune = true;
            }
            generationCts?.Cancel();
            IsGenerating = false;
            NotifyChanged();
        }

        public async Task GenerateNewPuzzleAsync()
        {
            if (IsGenerating) return;

            generationCts?.Cancel();
            generationCts = new CancellationTokenSource();
            var token = generationCts.Token;

            IsGenerating = true;
            GenerationStartTime = DateTime.UtcNow;
            IsSolved = false;
            PrunedProgress = 0;
            TotalIntersections = 0;
            NotifyChanged();

            // Allow UI to render the loading overlay before computation begins
            await Task.Delay(20);

            int width, height;
            if (!PuzzleHelper.ParseSize(SizeText, CurrentType, out width, out height))
            {
                PuzzleHelper.GetDefaultSize(CurrentType, out width, out height);
                SizeText = $"{width}x{height}";
            }

            try
            {
                Mesh mesh = PuzzleHelper.MakeMesh(width, height, CurrentType, Difficulty);
                generatingMesh = mesh;
                TotalIntersections = mesh.Intersections.Count;
                NotifyChanged();

                var progress = new Progress<int>(p =>
                {
                    PrunedProgress = p;
                    NotifyChanged();
                });

                await mesh.GenerateAsync(progress, token);

                if (!token.IsCancellationRequested && !mesh.AbortPrune)
                {
                    CurrentMesh = mesh;
                    MarkedEdges.Clear();
                    UndoTree = new UndoTree();
                    StartTime = DateTime.UtcNow;
                    IsSolved = false;
                }
            }
            catch (OperationCanceledException)
            {
                // Gracefully cancelled by user
            }
            finally
            {
                generatingMesh = null;
                generationCts = null;
                GenerationStartTime = null;
                IsGenerating = false;
                NotifyChanged();
            }
        }

        public void ToggleEdge(int edgeIndex, bool isAlternative)
        {
            if (CurrentMesh == null || edgeIndex < 0 || edgeIndex >= CurrentMesh.Edges.Count)
                return;

            if (MarkedEdges.Contains(edgeIndex))
                return;

            var current = CurrentMesh.Edges[edgeIndex].State;
            var next = !isAlternative
                ? (current == EdgeState.Empty ? EdgeState.Filled : (current == EdgeState.Filled ? EdgeState.Excluded : EdgeState.Empty))
                : (current == EdgeState.Empty ? EdgeState.Excluded : (current == EdgeState.Excluded ? EdgeState.Filled : EdgeState.Empty));

            var action = new PuzzleSetEdgeStateAction(CurrentMesh, edgeIndex, next, DisallowFalseMove);
            if (UndoTree.Do(action))
            {
                CheckSolution();
                NotifyChanged();
            }
        }

        public void SetEdgeDirect(int edgeIndex, EdgeState targetState)
        {
            if (CurrentMesh == null || edgeIndex < 0 || edgeIndex >= CurrentMesh.Edges.Count)
                return;

            if (MarkedEdges.Contains(edgeIndex))
                return;

            if (CurrentMesh.Edges[edgeIndex].State == targetState)
                return;

            var action = new PuzzleSetEdgeStateAction(CurrentMesh, edgeIndex, targetState, DisallowFalseMove);
            if (UndoTree.Do(action))
            {
                CheckSolution();
                NotifyChanged();
            }
        }

        public void ToggleCellColor(int cellIndex, bool isAlternative)
        {
            if (CurrentMesh == null || cellIndex < 0 || cellIndex >= CurrentMesh.Cells.Count)
                return;

            var action = new PuzzleCellColorAction(CurrentMesh, cellIndex, isAlternative);
            if (UndoTree.Do(action))
            {
                NotifyChanged();
            }
        }

        public void Undo()
        {
            if (UndoTree.CanUndo)
            {
                UndoTree.Undo();
                CheckSolution();
                NotifyChanged();
            }
        }

        public bool UndoAndForget()
        {
            if (UndoTree.CanUndo && UndoTree.UndoAndForget())
            {
                CheckSolution();
                NotifyChanged();
                return true;
            }
            return false;
        }

        public void Redo()
        {
            if (UndoTree.CanRedo)
            {
                UndoTree.Redo();
                CheckSolution();
                NotifyChanged();
            }
        }

        public void Fix()
        {
            if (CurrentMesh == null) return;
            MarkedEdges.Clear();
            for (int i = 0; i < CurrentMesh.Edges.Count; i++)
            {
                if (CurrentMesh.Edges[i].State != EdgeState.Empty)
                    MarkedEdges.Add(i);
            }
            UndoTree.Mark();
            NotifyChanged();
        }

        public void Unfix()
        {
            MarkedEdges.Clear();
            UndoTree.ClearMark();
            CheckSolution();
            NotifyChanged();
        }

        public void RevertToFix()
        {
            UndoTree.RevertToMark();
            CheckSolution();
            NotifyChanged();
        }

        public void CheckSolution()
        {
            if (CurrentMesh == null) return;
            bool currentlySolved = PuzzleHelper.CheckIsSolved(CurrentMesh);
            if (currentlySolved && !IsSolved)
            {
                IsSolved = true;
                Fix(); // Auto-fix solution upon victory
            }
            else if (!currentlySolved && IsSolved)
            {
                IsSolved = false;
            }
        }

        public SavedGameState? ExportState()
        {
            if (CurrentMesh == null) return null;

            int[] edgeStates = new int[CurrentMesh.Edges.Count];
            for (int i = 0; i < CurrentMesh.Edges.Count; i++)
            {
                edgeStates[i] = (int)CurrentMesh.Edges[i].State;
            }

            int[] cellColors = new int[CurrentMesh.Cells.Count];
            for (int i = 0; i < CurrentMesh.Cells.Count; i++)
            {
                cellColors[i] = CurrentMesh.Cells[i].Color;
            }

            PuzzleHelper.ParseSize(SizeText, CurrentType, out int w, out int h);
            string puzzleEncoded = PuzzleCodec.Encode(CurrentMesh, w, h);

            return new SavedGameState
            {
                Puzzle = puzzleEncoded,
                EdgeStates = edgeStates,
                CellColors = cellColors,
                MarkedEdges = MarkedEdges.ToArray(),
                ElapsedSeconds = ElapsedTime.TotalSeconds,
                IsSolved = IsSolved,
                DisallowFalseMove = DisallowFalseMove,
                ShowCellColors = ShowCellColors,
                Type = CurrentType.ToString(),
                Size = SizeText,
                Difficulty = Difficulty
            };
        }

        public bool RestoreState(SavedGameState state)
        {
            if (string.IsNullOrEmpty(state.Puzzle)) return false;

            Mesh? mesh = null;
            if (PuzzleCodec.TryDecode(state.Puzzle, out var decodedMesh, out int w, out int h, out var decodedType))
            {
                mesh = decodedMesh;
            }
            else
            {
                var lines = state.Puzzle.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                mesh = new Mesh(0, 0, MeshType.Square);
                if (!mesh.LoadFromText(lines)) return false;
            }

            if (mesh == null) return false;

            CurrentMesh = mesh;
            CurrentType = mesh.MeshType;
            if (!string.IsNullOrEmpty(state.Size)) SizeText = state.Size;
            Difficulty = state.Difficulty;
            DisallowFalseMove = state.DisallowFalseMove;
            ShowCellColors = state.ShowCellColors;
            TotalIntersections = mesh.Intersections.Count;
            PrunedProgress = 0;
            IsGenerating = false;
            UndoTree = new UndoTree();
            MarkedEdges.Clear();

            if (state.EdgeStates != null && state.EdgeStates.Length == mesh.Edges.Count)
            {
                for (int i = 0; i < state.EdgeStates.Length; i++)
                {
                    mesh.Edges[i].State = (EdgeState)state.EdgeStates[i];
                }
            }

            if (state.CellColors != null && state.CellColors.Length == mesh.Cells.Count)
            {
                for (int i = 0; i < state.CellColors.Length; i++)
                {
                    mesh.Cells[i].Color = state.CellColors[i];
                }
            }

            if (state.MarkedEdges != null)
            {
                foreach (int m in state.MarkedEdges)
                {
                    if (m >= 0 && m < mesh.Edges.Count)
                        MarkedEdges.Add(m);
                }
            }

            if (state.ElapsedSeconds > 0)
            {
                StartTime = DateTime.UtcNow.AddSeconds(-state.ElapsedSeconds);
            }
            else
            {
                StartTime = DateTime.UtcNow;
            }

            IsSolved = state.IsSolved;
            CheckSolution();
            NotifyChanged();
            return true;
        }

        public bool LoadPuzzleText(string puzzleText)
        {
            if (string.IsNullOrEmpty(puzzleText)) return false;

            Mesh? mesh = null;
            if (PuzzleCodec.TryDecode(puzzleText, out var decodedMesh, out int w, out int h, out var decodedType))
            {
                mesh = decodedMesh;
                SizeText = $"{w}x{h}";
            }
            else
            {
                var lines = puzzleText.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                mesh = new Mesh(0, 0, MeshType.Square);
                if (!mesh.LoadFromText(lines)) return false;
            }

            if (mesh == null) return false;

            CurrentMesh = mesh;
            CurrentType = mesh.MeshType;
            TotalIntersections = mesh.Intersections.Count;
            PrunedProgress = 0;
            IsGenerating = false;
            UndoTree = new UndoTree();
            MarkedEdges.Clear();
            StartTime = DateTime.UtcNow;
            IsSolved = false;

            CheckSolution();
            NotifyChanged();
            return true;
        }
    }

    public class SavedGameState
    {
        public string Puzzle { get; set; } = string.Empty;
        public int[]? EdgeStates { get; set; }
        public int[]? CellColors { get; set; }
        public int[]? MarkedEdges { get; set; }
        public double ElapsedSeconds { get; set; }
        public bool IsSolved { get; set; }
        public bool DisallowFalseMove { get; set; }
        public bool ShowCellColors { get; set; }
        public string Type { get; set; } = "Square";
        public string Size { get; set; } = "10x10";
        public int Difficulty { get; set; } = 1;
    }
}

