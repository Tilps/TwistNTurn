using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using TwistNTurn;

namespace TwistNTurn.Web.Services
{
    public sealed class GameService : IDisposable
    {
        private static readonly char[] LineSeparators = ['\r', '\n'];
        private const int MaxVisualSolveFrames = 300;

        public Mesh? CurrentMesh { get; private set; }
        public UndoTree UndoTree { get; private set; } = new();
        public HashSet<int> MarkedEdges { get; private set; } = new();

        public MeshType CurrentType { get; set; } = MeshType.Square;
        public string SizeText { get; set; } = "10x10";
        public int Difficulty { get; set; } = 1; // Easy

        public bool DisallowFalseMove { get; set; } = false;
        public bool ShowCellColors { get; set; } = false;

        private bool showEdgeColors = false;
        public bool ShowEdgeColors
        {
            get => showEdgeColors;
            set
            {
                showEdgeColors = value;
                if (!value) PendingGroupEdgeIndex = null;
            }
        }

        // Phase 2 Settings
        public int AutoMoveLevel { get; set; } = 0; // 0 = Off, 1..3 = Depth 1..3, 4 = Full Cascade
        public bool AutoMoveConsiderMultipleLoops { get; set; } = true;
        public bool AutoMoveUseColoring { get; set; } = false;
        public bool AutoStart { get; set; } = false;
        public MouseWheelMode WheelMode { get; set; } = MouseWheelMode.Zoom;
        public bool ShowClock { get; set; } = true;

        // Generator Overrides
        public double GenerateLengthFraction { get; set; } = 0.5;
        public double GenerateBoringFraction { get; set; } = 0.01;
        public GeneratorSolverMode SolverModeOverride { get; set; } = GeneratorSolverMode.Auto;

        public bool IsGenerating { get; private set; }
        public int PrunedProgress { get; private set; }
        public int TotalIntersections { get; private set; }

        public DateTime? StartTime { get; private set; }
        public TimeSpan ElapsedTime => StartTime.HasValue ? DateTime.UtcNow - StartTime.Value : TimeSpan.Zero;
        public bool IsSolved { get; private set; }
        public bool WasAutoSolved { get; private set; }

        // Phase 1 Visual Feedback & Hint State
        public int? InvalidEdgeIndex { get; private set; }
        public int? HintEdgeIndex { get; private set; }
        public EdgeState HintTargetState { get; private set; }
        public bool HintIsApplied { get; private set; }
        public string? HintMessage { get; private set; }

        // Phase 4 Loop Grouping & Inspector State
        public int? PendingGroupEdgeIndex { get; private set; }
        public PuzzleAnalysisResult? CurrentAnalysis { get; private set; }

        // Phase 3 Visual Solver State
        private readonly List<VisualSolveStep> visualSolveSteps = [];
        private Mesh? visualSolveSolutionMesh;
        public IReadOnlyList<VisualSolveStep> VisualSolveSteps => visualSolveSteps;
        public int VisualSolveTotalSteps => visualSolveSteps.Count;
        public int VisualSolveStepIndex { get; private set; }
        public bool IsVisualSolving => visualSolveSteps.Count > 0;
        public VisualSolveStep? CurrentVisualSolveStep => IsVisualSolving && VisualSolveStepIndex >= 0 && VisualSolveStepIndex < visualSolveSteps.Count
            ? visualSolveSteps[VisualSolveStepIndex]
            : null;

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

        private CancellationTokenSource? invalidFlashCts;
        private CancellationTokenSource? hintFlashCts;
        private CancellationTokenSource? hintMessageCts;

        private static CancellationToken ResetTokenSource(ref CancellationTokenSource? cts)
        {
            cts?.Cancel();
            cts?.Dispose();
            cts = new CancellationTokenSource();
            return cts.Token;
        }

        private async Task ClearAfterDelayAsync(int delayMs, CancellationToken token, Action clearAction)
        {
            try
            {
                await Task.Delay(delayMs, token);
                if (!token.IsCancellationRequested)
                {
                    clearAction();
                }
            }
            catch (OperationCanceledException)
            {
                // Expected when superseded or cleared early
            }
        }

        private void ClearTransientBoardState()
        {
            invalidFlashCts?.Cancel();
            hintFlashCts?.Cancel();
            hintMessageCts?.Cancel();
            InvalidEdgeIndex = null;
            HintEdgeIndex = null;
            HintMessage = null;
            PendingGroupEdgeIndex = null;
            CurrentAnalysis = null;
            preVictoryMarkedEdges = null;
            preVictoryUndoMark = null;
            preVictoryUndoTree = null;
            DismissVisualSolveInternal();
        }

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

            var token = ResetTokenSource(ref generationCts);

            IsGenerating = true;
            GenerationStartTime = DateTime.UtcNow;
            IsSolved = false;
            WasAutoSolved = false;
            ClearTransientBoardState();
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
                ApplyGeneratorOverrides(mesh);
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
                    ApplyAutoStartIfEnabled(mesh);
                    StartTime = DateTime.UtcNow;
                    IsSolved = false;
                    WasAutoSolved = false;
                    ClearTransientBoardState();
                    CheckSolution();
                }
            }
            catch (OperationCanceledException)
            {
                // Gracefully cancelled by user
            }
            finally
            {
                generatingMesh = null;
                generationCts?.Dispose();
                generationCts = null;
                GenerationStartTime = null;
                IsGenerating = false;
                NotifyChanged();
            }
        }

        private int GetEffectiveAutoMoveDepth(Mesh mesh)
        {
            return AutoMoveLevel switch
            {
                1 => 1,
                2 => 2,
                3 => 3,
                >= 4 => mesh.Edges.Count + 1,
                _ => 0
            };
        }

        public void ConfigureMeshForInteractiveMove(Mesh mesh)
        {
            mesh.ConsiderMultipleLoops = AutoMoveConsiderMultipleLoops;
            mesh.UseColoring = AutoMoveUseColoring;
            mesh.UseCellPairs = false;
            mesh.UseCellPairsTopLevel = false;
            mesh.UseEdgeRestricts = false;
            mesh.UseCellColoring = false;
            mesh.UseDerivedColoring = false;
            mesh.UseMerging = false;
            mesh.ColoringCheats = false;
        }

        public void ApplyGeneratorOverrides(Mesh mesh)
        {
            mesh.GenerateLengthFraction = GenerateLengthFraction;
            mesh.GenerateBoringFraction = GenerateBoringFraction;
            if (SolverModeOverride == GeneratorSolverMode.Iterative)
                mesh.SolverMethod = SolverMethod.Iterative;
            else if (SolverModeOverride == GeneratorSolverMode.Recursive)
                mesh.SolverMethod = SolverMethod.Recursive;
        }

        public void ApplyAutoStartIfEnabled(Mesh mesh)
        {
            if (!AutoStart) return;
            ConfigureMeshForInteractiveMove(mesh);
            mesh.PerformStart([]);
        }

        public UserPreferences ExportPreferences() => new()
        {
            DisallowFalseMove = DisallowFalseMove,
            ShowCellColors = ShowCellColors,
            ShowEdgeColors = ShowEdgeColors,
            AutoMoveLevel = AutoMoveLevel,
            AutoMoveConsiderMultipleLoops = AutoMoveConsiderMultipleLoops,
            AutoMoveUseColoring = AutoMoveUseColoring,
            AutoStart = AutoStart,
            WheelMode = WheelMode,
            ShowClock = ShowClock,
            GenerateLengthFraction = GenerateLengthFraction,
            GenerateBoringFraction = GenerateBoringFraction,
            SolverModeOverride = SolverModeOverride
        };

        public void ApplyPreferences(UserPreferences prefs)
        {
            DisallowFalseMove = prefs.DisallowFalseMove;
            ShowCellColors = prefs.ShowCellColors;
            ShowEdgeColors = prefs.ShowEdgeColors;
            AutoMoveLevel = prefs.AutoMoveLevel;
            AutoMoveConsiderMultipleLoops = prefs.AutoMoveConsiderMultipleLoops;
            AutoMoveUseColoring = prefs.AutoMoveUseColoring;
            AutoStart = prefs.AutoStart;
            WheelMode = Enum.IsDefined(prefs.WheelMode) ? prefs.WheelMode : MouseWheelMode.Zoom;
            ShowClock = prefs.ShowClock;
            GenerateLengthFraction = prefs.GenerateLengthFraction;
            GenerateBoringFraction = prefs.GenerateBoringFraction;
            SolverModeOverride = Enum.IsDefined(prefs.SolverModeOverride) ? prefs.SolverModeOverride : GeneratorSolverMode.Auto;
            NotifyChanged();
        }

        public void ResetPreferencesToDefaults()
        {
            ApplyPreferences(new UserPreferences());
        }

        public bool HandleEdgeGroupingClick(int edgeIndex, bool sameLoop)
        {
            if (CurrentMesh == null || IsGenerating || edgeIndex < 0 || edgeIndex >= CurrentMesh.Edges.Count)
                return false;

            bool dismissedVisual = DismissVisualSolveInternal();
            if (!ShowEdgeColors || MarkedEdges.Contains(edgeIndex))
            {
                if (dismissedVisual) NotifyChanged();
                return false;
            }

            if (PendingGroupEdgeIndex is not int firstEdge)
            {
                PendingGroupEdgeIndex = edgeIndex;
                NotifyChanged();
                return true;
            }

            PendingGroupEdgeIndex = null;
            if (firstEdge == edgeIndex)
            {
                NotifyChanged();
                return true;
            }

            var action = new ColorJoinAction(CurrentMesh, firstEdge, edgeIndex, sameLoop);
            if (UndoTree.Do(action))
            {
                InvalidEdgeIndex = null;
                NotifyChanged();
                return true;
            }

            FlashInvalidEdge(edgeIndex);
            NotifyChanged();
            return false;
        }

        public void ClearPendingEdgeGroup()
        {
            if (PendingGroupEdgeIndex.HasValue)
            {
                PendingGroupEdgeIndex = null;
                NotifyChanged();
            }
        }

        public bool CanExportSquareFormats => PuzzleFormatHelper.CanExportSquareFormats(CurrentMesh);

        public string ExportTwistText() => CurrentMesh != null ? PuzzleFormatHelper.ExportTwistText(CurrentMesh) : string.Empty;

        public string? ExportPearlDescriptor() => CurrentMesh != null ? PuzzleFormatHelper.ExportPearlDescriptor(CurrentMesh) : null;

        public string? ExportSquareAscii() => CurrentMesh != null ? PuzzleFormatHelper.ExportSquareAscii(CurrentMesh) : null;

        public bool ImportPuzzle(string rawInput, out string detectedFormat)
        {
            detectedFormat = string.Empty;
            if (IsGenerating || !PuzzleFormatHelper.TryParsePuzzle(rawInput, out Mesh? mesh, out detectedFormat, out int w, out int h) || mesh == null)
                return false;

            mesh.RecalculateCounts();
            CurrentMesh = mesh;
            CurrentType = mesh.MeshType;
            if (w > 0 && h > 0)
            {
                SizeText = $"{w}x{h}";
            }
            TotalIntersections = mesh.Intersections.Count;
            PrunedProgress = 0;
            IsGenerating = false;
            UndoTree = new UndoTree();
            MarkedEdges.Clear();
            StartTime = DateTime.UtcNow;
            IsSolved = false;
            WasAutoSolved = false;
            ClearTransientBoardState();

            CheckSolution();
            NotifyChanged();
            return true;
        }

        private static Mesh CreateFullSolverMesh(Mesh source)
        {
            var solver = new Mesh(source);
            solver.ConsiderMultipleLoops = true;
            solver.IterativeRecMaxDepth = 1;
            solver.SolverMethod = SolverMethod.Recursive;
            solver.UseCellColoring = true;
            solver.UseCellColoringTrials = true;
            solver.UseColoring = true;
            solver.UseEdgeRestricts = true;
            solver.UseDerivedColoring = true;
            solver.UseMerging = true;
            solver.UseCellPairs = true;
            solver.UseCellPairsTopLevel = true;
            return solver;
        }

        public PuzzleAnalysisResult? AnalyzeCurrentPuzzle()
        {
            if (CurrentMesh == null || IsGenerating) return null;

            int whitePearls = 0;
            int blackPearls = 0;
            foreach (var inter in CurrentMesh.Intersections)
            {
                if (inter.Type == IntersType.Twist) whitePearls++;
                else if (inter.Type == IntersType.Turn) blackPearls++;
            }
            int totalIntersections = CurrentMesh.Intersections.Count;
            int edgeCount = CurrentMesh.Edges.Count;

            Mesh solver = CreateFullSolverMesh(CurrentMesh);
            solver.Clear();

            SolveState fullResult = solver.TrySolve();
            if (fullResult == SolveState.MultipleSolutions)
            {
                CurrentAnalysis = new PuzzleAnalysisResult("Multiple Solutions", "Non-unique", 0, whitePearls, blackPearls, totalIntersections, edgeCount);
                NotifyChanged();
                return CurrentAnalysis;
            }
            if (fullResult != SolveState.Solved)
            {
                CurrentAnalysis = new PuzzleAnalysisResult("No Solution", "Unsolvable", 0, whitePearls, blackPearls, totalIntersections, edgeCount);
                NotifyChanged();
                return CurrentAnalysis;
            }

            // RecursiveTrySolve pushes a base frame before checking trials, so branch count is Length - 1.
            int branchDepth = Math.Max(0, (solver.DepthPatern?.Length ?? 0) - 1);
            string tier = $"Recursive Search (F, Depth {branchDepth})";
            int? lookaheadDepth = null;

            solver.Clear();
            solver.SetRatingCodeOptions("S");
            if (solver.TrySolve() == SolveState.Solved)
            {
                int min = -1;
                int max = solver.Edges.Count;
                while (min < max - 1)
                {
                    int mid = (min + max) / 2;
                    solver.Clear();
                    solver.IterativeSolverDepth = mid;
                    if (solver.TrySolve() == SolveState.Solved)
                        max = mid;
                    else
                        min = mid;
                }
                lookaheadDepth = max;
                tier = max == 0 ? "Direct Deductions (S0)" : $"Iterative Lookahead (S{max})";
            }
            else
            {
                (string code, string label)[] tiers =
                [
                    ("SC", "Edge Coloring (SC)"),
                    ("SC+", "Derived Coloring (SC+)"),
                    ("SCE", "Edge Restrictions (SCE)"),
                    ("SCOM", "Cell Coloring & Merging (SCOM)")
                ];
                foreach (var (code, label) in tiers)
                {
                    solver.Clear();
                    solver.SetRatingCodeOptions(code);
                    if (solver.TrySolve() == SolveState.Solved)
                    {
                        tier = label;
                        break;
                    }
                }
            }

            CurrentAnalysis = new PuzzleAnalysisResult("Unique Solution", tier, branchDepth, whitePearls, blackPearls, totalIntersections, edgeCount, lookaheadDepth);
            NotifyChanged();
            return CurrentAnalysis;
        }

        public void ToggleEdge(int edgeIndex, bool isAlternative)
        {
            if (CurrentMesh == null || edgeIndex < 0 || edgeIndex >= CurrentMesh.Edges.Count)
                return;

            PendingGroupEdgeIndex = null;
            bool dismissedVisual = DismissVisualSolveInternal();

            if (MarkedEdges.Contains(edgeIndex))
            {
                if (dismissedVisual) NotifyChanged();
                return;
            }

            EdgeState next = PuzzleEdgeAction.NextState(CurrentMesh.Edges[edgeIndex].State, isAlternative);
            int autoMove = GetEffectiveAutoMoveDepth(CurrentMesh);
            if (autoMove > 0)
                ConfigureMeshForInteractiveMove(CurrentMesh);

            var action = new PuzzleSetEdgeStateAction(CurrentMesh, edgeIndex, next, DisallowFalseMove, autoMove);
            if (UndoTree.Do(action))
            {
                InvalidEdgeIndex = null;
                HintEdgeIndex = null;
                CheckSolution();
                NotifyChanged();
            }
            else
            {
                if (dismissedVisual) NotifyChanged();
                if (DisallowFalseMove)
                    FlashInvalidEdge(edgeIndex);
            }
        }

        public void SetEdgeDirect(int edgeIndex, EdgeState targetState)
        {
            if (CurrentMesh == null || edgeIndex < 0 || edgeIndex >= CurrentMesh.Edges.Count)
                return;

            bool dismissedVisual = DismissVisualSolveInternal();

            if (MarkedEdges.Contains(edgeIndex) || CurrentMesh.Edges[edgeIndex].State == targetState)
            {
                if (dismissedVisual) NotifyChanged();
                return;
            }

            int autoMove = GetEffectiveAutoMoveDepth(CurrentMesh);
            if (autoMove > 0)
                ConfigureMeshForInteractiveMove(CurrentMesh);

            var action = new PuzzleSetEdgeStateAction(CurrentMesh, edgeIndex, targetState, DisallowFalseMove, autoMove);
            if (UndoTree.Do(action))
            {
                InvalidEdgeIndex = null;
                HintEdgeIndex = null;
                CheckSolution();
                NotifyChanged();
            }
            else
            {
                if (dismissedVisual) NotifyChanged();
                if (DisallowFalseMove)
                    FlashInvalidEdge(edgeIndex);
            }
        }

        public void ToggleCellColor(int cellIndex, bool isAlternative)
        {
            if (CurrentMesh == null || cellIndex < 0 || cellIndex >= CurrentMesh.Cells.Count)
                return;

            bool dismissedVisual = DismissVisualSolveInternal();
            var action = new PuzzleCellColorAction(CurrentMesh, cellIndex, isAlternative);
            if (UndoTree.Do(action))
            {
                NotifyChanged();
            }
            else if (dismissedVisual)
            {
                NotifyChanged();
            }
        }

        public void Undo()
        {
            bool dismissedVisual = DismissVisualSolveInternal();
            if (UndoTree.CanUndo)
            {
                UndoTree.Undo();
                WasAutoSolved = false;
                CheckSolution();
                NotifyChanged();
            }
            else if (dismissedVisual)
            {
                NotifyChanged();
            }
        }

        public bool UndoAndForget()
        {
            bool dismissedVisual = DismissVisualSolveInternal();
            if (UndoTree.CanUndo && UndoTree.UndoAndForget())
            {
                WasAutoSolved = false;
                CheckSolution();
                NotifyChanged();
                return true;
            }
            if (dismissedVisual) NotifyChanged();
            return false;
        }

        public void Redo()
        {
            bool dismissedVisual = DismissVisualSolveInternal();
            if (UndoTree.CanRedo)
            {
                UndoTree.Redo();
                WasAutoSolved = UndoTree.CurrentAction is PuzzleBatchSolveAction;
                CheckSolution();
                NotifyChanged();
            }
            else if (dismissedVisual)
            {
                NotifyChanged();
            }
        }

        public void Fix()
        {
            if (CurrentMesh == null) return;
            DismissVisualSolveInternal();
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
            DismissVisualSolveInternal();
            MarkedEdges.Clear();
            UndoTree.ClearMark();
            CheckSolution();
            NotifyChanged();
        }

        public void RevertToFix()
        {
            DismissVisualSolveInternal();
            UndoTree.RevertToMark();
            CheckSolution();
            NotifyChanged();
        }

        // Checkpoint state captured just before the victory auto-Fix, so it can be
        // restored if the puzzle becomes unsolved again (e.g. via Undo).
        private int[]? preVictoryMarkedEdges;
        private object? preVictoryUndoMark;
        private UndoTree? preVictoryUndoTree;

        public void CheckSolution()
        {
            if (CurrentMesh == null) return;
            bool currentlySolved = PuzzleHelper.CheckIsSolved(CurrentMesh);
            if (currentlySolved && !IsSolved)
            {
                preVictoryMarkedEdges = MarkedEdges.ToArray();
                preVictoryUndoMark = UndoTree.GetMarkedDirect();
                preVictoryUndoTree = UndoTree;
                IsSolved = true;
                Fix(); // Auto-fix solution upon victory
            }
            else if (!currentlySolved && IsSolved)
            {
                IsSolved = false;
                if (preVictoryMarkedEdges != null && ReferenceEquals(preVictoryUndoTree, UndoTree))
                {
                    MarkedEdges = [.. preVictoryMarkedEdges];
                    UndoTree.SetMarkedDirect(preVictoryUndoMark);
                    preVictoryMarkedEdges = null;
                    preVictoryUndoMark = null;
                    preVictoryUndoTree = null;
                }
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
                WasAutoSolved = WasAutoSolved,
                DisallowFalseMove = DisallowFalseMove,
                ShowCellColors = ShowCellColors,
                Type = CurrentType.ToString(),
                Size = SizeText,
                Difficulty = Difficulty
            };
        }

        private static bool TryDecodeOrLoadMesh(string data, out Mesh? mesh, out bool isDecoded, out int width, out int height)
        {
            if (PuzzleCodec.TryDecode(data, out mesh, out width, out height, out _))
            {
                isDecoded = true;
                return mesh != null;
            }

            isDecoded = false;
            var lines = data.Split(LineSeparators, StringSplitOptions.RemoveEmptyEntries);
            mesh = new Mesh(0, 0, MeshType.Square);
            if (!mesh.LoadFromText(lines))
            {
                mesh = null;
                return false;
            }
            return true;
        }

        public bool RestoreState(SavedGameState state)
        {
            if (string.IsNullOrEmpty(state.Puzzle)) return false;

            if (!TryDecodeOrLoadMesh(state.Puzzle, out var mesh, out _, out _, out _) || mesh == null)
                return false;

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
            mesh.RecalculateCounts();

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

            StartTime = state.ElapsedSeconds > 0
                ? DateTime.UtcNow.AddSeconds(-state.ElapsedSeconds)
                : DateTime.UtcNow;

            IsSolved = state.IsSolved;
            WasAutoSolved = state.IsSolved && state.WasAutoSolved;
            ClearTransientBoardState();
            CheckSolution();
            NotifyChanged();
            return true;
        }

        public bool LoadPuzzleText(string puzzleText)
        {
            if (string.IsNullOrEmpty(puzzleText)) return false;

            if (!TryDecodeOrLoadMesh(puzzleText, out var mesh, out bool isDecoded, out int w, out int h) || mesh == null)
                return false;

            if (isDecoded && w > 0 && h > 0)
            {
                SizeText = $"{w}x{h}";
            }

            mesh.RecalculateCounts();
            CurrentMesh = mesh;
            CurrentType = mesh.MeshType;
            TotalIntersections = mesh.Intersections.Count;
            PrunedProgress = 0;
            IsGenerating = false;
            UndoTree = new UndoTree();
            MarkedEdges.Clear();
            StartTime = DateTime.UtcNow;
            IsSolved = false;
            WasAutoSolved = false;
            ClearTransientBoardState();

            CheckSolution();
            NotifyChanged();
            return true;
        }

        public void FlashInvalidEdge(int edgeIndex)
        {
            var token = ResetTokenSource(ref invalidFlashCts);

            InvalidEdgeIndex = edgeIndex;
            NotifyChanged();

            _ = ClearAfterDelayAsync(600, token, () =>
            {
                if (InvalidEdgeIndex == edgeIndex)
                {
                    InvalidEdgeIndex = null;
                    NotifyChanged();
                }
            });
        }

        public void ClearHint()
        {
            hintFlashCts?.Cancel();
            hintMessageCts?.Cancel();
            HintEdgeIndex = null;
            HintMessage = null;
            NotifyChanged();
        }

        public void ResetCurrentPuzzle()
        {
            if (CurrentMesh == null || IsGenerating) return;

            CurrentMesh.Clear();
            MarkedEdges.Clear();
            UndoTree = new UndoTree();
            ApplyAutoStartIfEnabled(CurrentMesh);
            StartTime = DateTime.UtcNow;
            IsSolved = false;
            WasAutoSolved = false;
            ClearTransientBoardState();
            CheckSolution();
            NotifyChanged();
        }

        private static EdgeState[] CaptureEdgeStates(Mesh mesh)
        {
            var states = new EdgeState[mesh.Edges.Count];
            for (int i = 0; i < mesh.Edges.Count; i++)
            {
                states[i] = mesh.Edges[i].State;
            }
            return states;
        }

        private bool DismissVisualSolveInternal()
        {
            if (visualSolveSteps.Count == 0 && visualSolveSolutionMesh == null)
                return false;

            visualSolveSteps.Clear();
            visualSolveSolutionMesh = null;
            VisualSolveStepIndex = 0;
            return true;
        }

        public void DismissVisualSolve()
        {
            if (DismissVisualSolveInternal())
                NotifyChanged();
        }

        public bool StartVisualSolve()
        {
            if (CurrentMesh == null || IsGenerating) return false;

            DismissVisualSolveInternal();
            ClearHint();

            Mesh solver = CreateFullSolverMesh(CurrentMesh);
            solver.Clear();
            solver.SuperSlowMo = true;

            visualSolveSteps.Add(new VisualSolveStep(CaptureEdgeStates(solver), "Starting solver..."));

            void OnUpdate(object sender, MeshChangeUpdateEventArgs args)
            {
                if (visualSolveSteps.Count >= MaxVisualSolveFrames - 1)
                    return;

                EdgeState[] snapshot = CaptureEdgeStates(args.CurrentMesh);
                if (snapshot.AsSpan().SequenceEqual(visualSolveSteps[^1].EdgeStates))
                    return;

                string status;
                bool isBacktrack = false;
                bool isBranching = false;

                if (args.Starting)
                {
                    status = "Testing branch...";
                    isBranching = true;
                }
                else if (!args.SuccessfulAttempt && args.JustDone != null)
                {
                    status = "Contradiction found — backtracking";
                    isBacktrack = true;
                }
                else if (args.JustDone != null)
                {
                    status = "Branch consistent — evaluating opposite";
                    isBranching = true;
                }
                else if (!args.SuccessfulAttempt)
                {
                    status = "Propagating deduction wave...";
                }
                else
                {
                    status = "Applying deduced moves";
                }

                visualSolveSteps.Add(new VisualSolveStep(snapshot, status, isBacktrack, isBranching));
            }

            solver.MeshChangeUpdate += OnUpdate;
            SolveState res;
            try
            {
                res = solver.TrySolve();
            }
            finally
            {
                solver.MeshChangeUpdate -= OnUpdate;
            }

            if ((res == SolveState.Solved || res == SolveState.MultipleSolutions) && solver.SolutionFound != null)
            {
                EdgeState[] finalStates = CaptureEdgeStates(solver.SolutionFound);
                visualSolveSolutionMesh = solver.SolutionFound;

                if (!finalStates.AsSpan().SequenceEqual(visualSolveSteps[^1].EdgeStates))
                {
                    visualSolveSteps.Add(new VisualSolveStep(finalStates, "Solution complete!"));
                }
                else
                {
                    visualSolveSteps[^1] = visualSolveSteps[^1] with { Status = "Solution complete!", IsBacktrack = false, IsBranching = false };
                }

                VisualSolveStepIndex = 0;
                NotifyChanged();
                return true;
            }

            DismissVisualSolveInternal();
            NotifyChanged();
            return false;
        }

        public void SetVisualSolveStep(int index)
        {
            if (!IsVisualSolving) return;
            int clamped = Math.Clamp(index, 0, visualSolveSteps.Count - 1);
            if (clamped == VisualSolveStepIndex) return;
            VisualSolveStepIndex = clamped;
            NotifyChanged();
        }

        public bool StepVisualSolveForward()
        {
            if (!IsVisualSolving || VisualSolveStepIndex >= visualSolveSteps.Count - 1)
                return false;
            VisualSolveStepIndex++;
            NotifyChanged();
            return true;
        }

        public bool StepVisualSolveBackward()
        {
            if (!IsVisualSolving || VisualSolveStepIndex <= 0)
                return false;
            VisualSolveStepIndex--;
            NotifyChanged();
            return true;
        }

        public bool ApplyVisualSolveSolution()
        {
            if (!IsVisualSolving || CurrentMesh == null || visualSolveSolutionMesh == null)
                return false;

            var action = new PuzzleBatchSolveAction(CurrentMesh, visualSolveSolutionMesh);
            if (UndoTree.Do(action))
            {
                WasAutoSolved = true;
                DismissVisualSolveInternal();
                CheckSolution();
                NotifyChanged();
                return true;
            }

            return false;
        }

        public bool SolvePuzzle()
        {
            if (CurrentMesh == null || IsGenerating) return false;

            DismissVisualSolveInternal();
            Mesh solver = CreateFullSolverMesh(CurrentMesh);
            solver.Clear();

            SolveState res = solver.TrySolve();
            if ((res == SolveState.Solved || res == SolveState.MultipleSolutions) && solver.SolutionFound != null)
            {
                var action = new PuzzleBatchSolveAction(CurrentMesh, solver.SolutionFound);
                if (UndoTree.Do(action))
                {
                    WasAutoSolved = true;
                    CheckSolution();
                    NotifyChanged();
                    return true;
                }
            }
            return false;
        }

        public bool RequestHint(bool applyMove)
        {
            if (CurrentMesh == null || IsGenerating) return false;
            DismissVisualSolveInternal();

            var (edgeIdx, targetState, message) = ComputeHint();
            HintMessage = message;

            if (edgeIdx.HasValue)
            {
                HintEdgeIndex = edgeIdx.Value;
                HintTargetState = targetState;
                HintIsApplied = applyMove;
                if (applyMove)
                {
                    var action = new PuzzleSetEdgeStateAction(CurrentMesh, edgeIdx.Value, targetState, DisallowFalseMove);
                    UndoTree.Do(action);
                    CheckSolution();
                }
                else
                {
                    // Don't give away whether it's a line or a cross.
                    HintMessage = "The highlighted edge can be determined from what's around it.";
                }

                var flashToken = ResetTokenSource(ref hintFlashCts);
                int delayMs = applyMove ? 1500 : 3500;
                int capturedEdge = edgeIdx.Value;
                _ = ClearAfterDelayAsync(delayMs, flashToken, () =>
                {
                    if (HintEdgeIndex == capturedEdge)
                    {
                        HintEdgeIndex = null;
                        NotifyChanged();
                    }
                });
            }

            var msgToken = ResetTokenSource(ref hintMessageCts);
            _ = ClearAfterDelayAsync(5000, msgToken, () =>
            {
                HintMessage = null;
                NotifyChanged();
            });

            NotifyChanged();
            return edgeIdx.HasValue;
        }

        public (int? EdgeIndex, EdgeState TargetState, string Message) ComputeHint()
        {
            if (CurrentMesh == null)
                return (null, EdgeState.Empty, "No puzzle active.");

            if (IsSolved)
                return (null, EdgeState.Empty, "The puzzle is already solved!");

            // 1. Check for immediate Masyu contradictions on the board
            CurrentMesh.RecalculateCounts();
            for (int v = 0; v < CurrentMesh.Intersections.Count; v++)
            {
                var inter = CurrentMesh.Intersections[v];
                if (inter.FilledCount > 2)
                {
                    return (null, EdgeState.Empty, "Contradiction: An intersection has more than 2 connected lines.");
                }
                if (inter.FilledCount == 1 && inter.Edges.Count - inter.ExcludedCount < 2)
                {
                    return (null, EdgeState.Empty, "Contradiction: A loop end at an intersection has no available exit path.");
                }
                if (inter.Type != IntersType.Unknown && inter.Edges.Count - inter.ExcludedCount < 2)
                {
                    return (null, EdgeState.Empty, "Contradiction: A pearl intersection does not have enough available edges to pass through.");
                }
                if (inter.Type != IntersType.Unknown && inter.FilledCount == 2 && TryCheckLocalStraight(CurrentMesh, inter, out bool isStraight))
                {
                    if (inter.Type == IntersType.Twist && !isStraight)
                    {
                        return (null, EdgeState.Empty, "Contradiction: Loop cannot turn at a white pearl.");
                    }
                    if (inter.Type == IntersType.Turn && isStraight)
                    {
                        return (null, EdgeState.Empty, "Contradiction: Loop cannot go straight through a black pearl.");
                    }
                }
            }

            // 2. Clone mesh and run iterative deduction pass
            Mesh hintMesh = new Mesh(CurrentMesh);
            hintMesh.ConsiderMultipleLoops = true;
            hintMesh.IterativeSolverDepth = 0;
            hintMesh.IterativeRecMaxDepth = 1;
            hintMesh.UseColoring = false;
            hintMesh.UseCellPairs = false;
            hintMesh.UseCellPairsTopLevel = false;
            hintMesh.UseEdgeRestricts = false;
            hintMesh.UseCellColoring = false;
            hintMesh.UseDerivedColoring = false;
            hintMesh.UseMerging = false;

            List<IAction> changes = [];
            bool startOk = hintMesh.PerformStart(changes);
            if (!startOk)
            {
                return (null, EdgeState.Empty, "Contradiction detected in current board state.");
            }

            SetAction? action = changes.OfType<SetAction>().FirstOrDefault(a =>
                a.EdgeIndex >= 0 && a.EdgeIndex < CurrentMesh.Edges.Count &&
                CurrentMesh.Edges[a.EdgeIndex].State != a.EdgeState);

            while (action == null && hintMesh.IterativeSolverDepth <= 10)
            {
                changes.Clear();
                hintMesh.GetSomething(changes);
                action = changes.OfType<SetAction>().FirstOrDefault(a =>
                    a.EdgeIndex >= 0 && a.EdgeIndex < CurrentMesh.Edges.Count &&
                    CurrentMesh.Edges[a.EdgeIndex].State != a.EdgeState);
                if (action != null)
                    break;
                hintMesh.IterativeSolverDepth++;
            }

            if (action != null)
            {
                int edgeIdx = action.EdgeIndex;
                EdgeState state = action.EdgeState;
                string reason = ExplainDeduction(CurrentMesh, edgeIdx, state);
                return (edgeIdx, state, reason);
            }

            // 3. Fallback: full solver search on current state
            Mesh fullSolverMesh = CreateFullSolverMesh(CurrentMesh);
            SolveState solveRes = fullSolverMesh.TrySolve();
            if ((solveRes == SolveState.Solved || solveRes == SolveState.MultipleSolutions) && fullSolverMesh.SolutionFound != null)
            {
                var sol = fullSolverMesh.SolutionFound;
                for (int i = 0; i < CurrentMesh.Edges.Count; i++)
                {
                    if (CurrentMesh.Edges[i].State == EdgeState.Empty && sol.Edges[i].State != EdgeState.Empty)
                    {
                        string stateStr = sol.Edges[i].State == EdgeState.Filled ? "line" : "cross (✕)";
                        return (i, sol.Edges[i].State, $"Advanced deduction: This edge must be a {stateStr}.");
                    }
                }
            }
            else if (solveRes == SolveState.NoSolutions)
            {
                return (null, EdgeState.Empty, "Contradiction: No valid solution exists from this board state.");
            }

            return (null, EdgeState.Empty, "No logical deduction found at this step.");
        }

        private static bool TryCheckLocalStraight(Mesh mesh, Intersection inter, out bool isStraight)
        {
            isStraight = false;
            int firstIdx = -1;
            int secondIdx = -1;
            for (int i = 0; i < inter.Edges.Count; i++)
            {
                if (mesh.Edges[inter.Edges[i]].State == EdgeState.Filled)
                {
                    if (firstIdx == -1) firstIdx = i;
                    else if (secondIdx == -1) secondIdx = i;
                    else return false;
                }
            }
            if (firstIdx < 0 || secondIdx < 0) return false;
            isStraight = inter.Opposite[firstIdx] == inter.Edges[secondIdx];
            return true;
        }

        private static string ExplainDeduction(Mesh mesh, int edgeIndex, EdgeState targetState)
        {
            Edge edge = mesh.Edges[edgeIndex];
            if (edge.Intersections != null)
            {
                foreach (int vIdx in edge.Intersections)
                {
                    if (vIdx < 0 || vIdx >= mesh.Intersections.Count) continue;
                    var inter = mesh.Intersections[vIdx];
                    if (inter.FilledCount == 2 && targetState == EdgeState.Excluded)
                    {
                        return "Intersection already has 2 connections: remaining edges must be excluded.";
                    }
                    if (inter.FilledCount == 1 && inter.Edges.Count - inter.ExcludedCount == 2 && targetState == EdgeState.Filled)
                    {
                        return "Loop entering intersection has only one exit path: must continue here.";
                    }
                    if (inter.FilledCount == 0 && inter.Edges.Count - inter.ExcludedCount < 2 && targetState == EdgeState.Excluded)
                    {
                        return "Dead end: loop cannot enter intersection without an exit path.";
                    }
                    if (inter.Type == IntersType.Twist)
                    {
                        return "White pearl rule: loop must pass straight through and turn immediately before or after.";
                    }
                    if (inter.Type == IntersType.Turn)
                    {
                        return "Black pearl rule: loop must turn at the pearl and extend straight on both sides.";
                    }
                }
            }

            return targetState == EdgeState.Filled
                ? "Logical deduction: This edge must be part of the loop."
                : "Logical deduction: This edge cannot be part of the loop.";
        }

        public void Dispose()
        {
            generationCts?.Cancel();
            generationCts?.Dispose();
            invalidFlashCts?.Cancel();
            invalidFlashCts?.Dispose();
            hintFlashCts?.Cancel();
            hintFlashCts?.Dispose();
            hintMessageCts?.Cancel();
            hintMessageCts?.Dispose();
        }
    }

    public enum MouseWheelMode
    {
        Zoom = 0,
        UndoRedo = 1
    }

    public enum GeneratorSolverMode
    {
        Auto = 0,
        Iterative = 1,
        Recursive = 2
    }

    public sealed record VisualSolveStep(
        EdgeState[] EdgeStates,
        string Status,
        bool IsBacktrack = false,
        bool IsBranching = false);

    public sealed record PuzzleAnalysisResult(
        string SolveStatus,
        string TechniqueTier,
        int BranchDepth,
        int WhitePearls,
        int BlackPearls,
        int TotalIntersections,
        int EdgeCount,
        int? LookaheadDepth = null)
    {
        public int ClueCount => WhitePearls + BlackPearls;
    }

    public class UserPreferences
    {
        public bool DisallowFalseMove { get; set; } = false;
        public bool ShowCellColors { get; set; } = false;
        public bool ShowEdgeColors { get; set; } = false;
        public int AutoMoveLevel { get; set; } = 0;
        public bool AutoMoveConsiderMultipleLoops { get; set; } = true;
        public bool AutoMoveUseColoring { get; set; } = false;
        public bool AutoStart { get; set; } = false;
        public MouseWheelMode WheelMode { get; set; } = MouseWheelMode.Zoom;
        public bool ShowClock { get; set; } = true;
        public double GenerateLengthFraction { get; set; } = 0.5;
        public double GenerateBoringFraction { get; set; } = 0.01;
        public GeneratorSolverMode SolverModeOverride { get; set; } = GeneratorSolverMode.Auto;
    }

    public class SavedGameState
    {
        public string Puzzle { get; set; } = string.Empty;
        public int[]? EdgeStates { get; set; }
        public int[]? CellColors { get; set; }
        public int[]? MarkedEdges { get; set; }
        public double ElapsedSeconds { get; set; }
        public bool IsSolved { get; set; }
        public bool WasAutoSolved { get; set; }
        public bool DisallowFalseMove { get; set; }
        public bool ShowCellColors { get; set; }
        public string Type { get; set; } = "Square";
        public string Size { get; set; } = "10x10";
        public int Difficulty { get; set; } = 1;
    }
}
