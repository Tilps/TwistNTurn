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

        // Active puzzle geometry metadata (independent of pending form edits)
        public int ActiveWidth { get; private set; }
        public int ActiveHeight { get; private set; }
        public MeshType ActiveType => CurrentMesh?.MeshType ?? CurrentType;

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
        private int autoMoveLevel = 0;
        public int AutoMoveLevel
        {
            get => autoMoveLevel;
            set => autoMoveLevel = Math.Clamp(value, 0, 4);
        }

        public bool AutoMoveConsiderMultipleLoops { get; set; } = true;
        public bool AutoMoveUseColoring { get; set; } = false;
        public bool AutoStart { get; set; } = false;
        public MouseWheelMode WheelMode { get; set; } = MouseWheelMode.Zoom;
        public bool ShowClock { get; set; } = true;

        // Generator Overrides
        private double generateLengthFraction = 0.5;
        public double GenerateLengthFraction
        {
            get => generateLengthFraction;
            set => generateLengthFraction = double.IsFinite(value) ? Math.Clamp(value, 0.25, 0.75) : 0.5;
        }

        private double generateBoringFraction = 0.01;
        public double GenerateBoringFraction
        {
            get => generateBoringFraction;
            set => generateBoringFraction = double.IsFinite(value) ? Math.Clamp(value, 0.0, 0.20) : 0.01;
        }

        public GeneratorSolverMode SolverModeOverride { get; set; } = GeneratorSolverMode.Auto;

        public bool IsGenerating { get; private set; }
        public bool IsBusy { get; private set; }
        public string? BusyMessage { get; private set; }
        public int PrunedProgress { get; private set; }
        public int TotalIntersections { get; private set; }

        private TimeSpan? solvedElapsed;
        public DateTime? StartTime { get; private set; }
        public TimeSpan ElapsedTime => solvedElapsed ?? (StartTime.HasValue ? DateTime.UtcNow - StartTime.Value : TimeSpan.Zero);
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
        public event Action? OnBoardPersistedStateChanged;

        public bool CanUndo => UndoTree.CanUndo;
        public bool CanRedo => UndoTree.CanRedo;

        public void NotifyChanged()
        {
            OnStateChanged?.Invoke();
            OnBoardPersistedStateChanged?.Invoke();
        }

        public void NotifyVisualChanged()
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
            var previous = cts;
            var next = new CancellationTokenSource();
            cts = next;
            if (previous != null)
            {
                previous.Cancel();
                previous.Dispose();
            }
            return next.Token;
        }

        private static void CancelAndDispose(ref CancellationTokenSource? cts)
        {
            var previous = cts;
            cts = null;
            if (previous != null)
            {
                previous.Cancel();
                previous.Dispose();
            }
        }

        private static async Task ClearAfterDelayAsync(int delayMs, Action clearAction, CancellationToken cancellationToken)
        {
            try
            {
                await Task.Delay(delayMs, cancellationToken);
                if (!cancellationToken.IsCancellationRequested)
                {
                    clearAction();
                }
            }
            catch (OperationCanceledException)
            {
                // Expected when superseded or cleared early
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Transient clear error: {ex.Message}");
            }
        }

        private void ClearHintInternal()
        {
            hintFlashCts?.Cancel();
            hintMessageCts?.Cancel();
            HintEdgeIndex = null;
            HintMessage = null;
        }

        private void ClearTransientBoardState()
        {
            invalidFlashCts?.Cancel();
            ClearHintInternal();
            InvalidEdgeIndex = null;
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
            GenerationStartTime = null;
            NotifyVisualChanged();
        }

        public async Task GenerateNewPuzzleAsync()
        {
            if (IsGenerating || IsBusy) return;

            var myCts = new CancellationTokenSource();
            var previous = generationCts;
            generationCts = myCts;
            if (previous != null)
            {
                previous.Cancel();
            }
            var token = myCts.Token;

            DismissVisualSolveInternal();
            IsGenerating = true;
            GenerationStartTime = DateTime.UtcNow;
            PrunedProgress = 0;
            TotalIntersections = 0;
            NotifyVisualChanged();

            bool generatedNewBoard = false;
            try
            {
                // Allow UI to render the loading overlay before computation begins
                await Task.Delay(20, token);

                MeshType targetType = CurrentType;
                int width, height;
                if (!PuzzleHelper.ParseSize(SizeText, targetType, out width, out height))
                {
                    PuzzleHelper.GetDefaultSize(targetType, out width, out height);
                    SizeText = $"{width}x{height}";
                }

                Mesh mesh = PuzzleHelper.MakeMesh(width, height, targetType, Difficulty);
                ApplyGeneratorOverrides(mesh);
                generatingMesh = mesh;
                TotalIntersections = mesh.Intersections.Count;
                NotifyVisualChanged();

                var progress = new Progress<int>(p =>
                {
                    if (ReferenceEquals(generationCts, myCts))
                    {
                        PrunedProgress = p;
                        NotifyVisualChanged();
                    }
                });

                await mesh.GenerateAsync(progress, token);

                if (!token.IsCancellationRequested && !mesh.AbortPrune && ReferenceEquals(generationCts, myCts))
                {
                    CurrentMesh = mesh;
                    ActiveWidth = width;
                    ActiveHeight = height;
                    MarkedEdges.Clear();
                    UndoTree = new UndoTree();
                    ApplyAutoStartIfEnabled(mesh);
                    StartTime = DateTime.UtcNow;
                    solvedElapsed = null;
                    IsSolved = false;
                    WasAutoSolved = false;
                    ClearTransientBoardState();
                    CheckSolution();
                    generatedNewBoard = true;
                }
            }
            catch (OperationCanceledException)
            {
                // Gracefully cancelled by user
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Puzzle generation error: {ex.Message}");
            }
            finally
            {
                if (ReferenceEquals(generationCts, myCts))
                {
                    generatingMesh = null;
                    generationCts = null;
                    GenerationStartTime = null;
                    IsGenerating = false;
                }
                myCts.Dispose();

                if (generatedNewBoard)
                    NotifyChanged();
                else
                    NotifyVisualChanged();
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
            NotifyVisualChanged();
        }

        public void ResetPreferencesToDefaults()
        {
            ApplyPreferences(new UserPreferences());
        }

        public bool HandleEdgeGroupingClick(int edgeIndex, bool sameLoop)
        {
            if (CurrentMesh == null || IsGenerating || IsBusy || edgeIndex < 0 || edgeIndex >= CurrentMesh.Edges.Count)
                return false;

            bool dismissedVisual = DismissVisualSolveInternal();
            if (!ShowEdgeColors || MarkedEdges.Contains(edgeIndex))
            {
                if (dismissedVisual) NotifyVisualChanged();
                return false;
            }

            if (PendingGroupEdgeIndex is not int firstEdge)
            {
                PendingGroupEdgeIndex = edgeIndex;
                NotifyVisualChanged();
                return true;
            }

            PendingGroupEdgeIndex = null;
            if (firstEdge == edgeIndex)
            {
                NotifyVisualChanged();
                return true;
            }

            var action = new ColorJoinAction(CurrentMesh, firstEdge, edgeIndex, sameLoop);
            if (UndoTree.Do(action))
            {
                InvalidEdgeIndex = null;
                NotifyVisualChanged();
                return true;
            }

            FlashInvalidEdge(edgeIndex);
            return false;
        }

        public void ClearPendingEdgeGroup()
        {
            if (PendingGroupEdgeIndex.HasValue)
            {
                PendingGroupEdgeIndex = null;
                NotifyVisualChanged();
            }
        }

        public bool CanExportSquareFormats => PuzzleFormatHelper.CanExportSquareFormats(CurrentMesh);

        public string ExportTwistText() => CurrentMesh != null ? PuzzleFormatHelper.ExportTwistText(CurrentMesh) : string.Empty;

        public string? ExportPearlDescriptor() => CurrentMesh != null ? PuzzleFormatHelper.ExportPearlDescriptor(CurrentMesh) : null;

        public string? ExportSquareAscii() => CurrentMesh != null ? PuzzleFormatHelper.ExportSquareAscii(CurrentMesh) : null;

        public string GetShareCode()
        {
            if (CurrentMesh == null) return string.Empty;

            if (ActiveWidth > 0 && ActiveHeight > 0)
            {
                return PuzzleCodec.Encode(CurrentMesh, ActiveWidth, ActiveHeight);
            }

            if (PuzzleFormatHelper.TryGetSquareDimensions(CurrentMesh, out int sw, out int sh))
            {
                return PuzzleCodec.Encode(CurrentMesh, sw, sh);
            }

            return PuzzleCompression.CompressToUrlSafe(CurrentMesh.SaveToString());
        }

        public bool ImportPuzzle(string rawInput, out string detectedFormat)
        {
            detectedFormat = string.Empty;
            if (IsGenerating || IsBusy || !PuzzleFormatHelper.TryParsePuzzle(rawInput, out Mesh? mesh, out detectedFormat, out int w, out int h) || mesh == null)
                return false;

            mesh.RecalculateCounts();
            CurrentMesh = mesh;
            CurrentType = mesh.MeshType;
            ActiveWidth = w;
            ActiveHeight = h;
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
            solvedElapsed = null;
            IsSolved = false;
            WasAutoSolved = false;
            ClearTransientBoardState();

            CheckSolution();
            NotifyChanged();
            return true;
        }

        public PuzzleAnalysisResult? AnalyzeCurrentPuzzle()
        {
            if (CurrentMesh == null || IsGenerating) return null;

            CurrentAnalysis = PuzzleHelper.AnalyzePuzzle(CurrentMesh);
            NotifyVisualChanged();
            return CurrentAnalysis;
        }

        public async Task<PuzzleAnalysisResult?> AnalyzeCurrentPuzzleAsync()
        {
            if (CurrentMesh == null || IsGenerating || IsBusy) return null;
            IsBusy = true;
            BusyMessage = "Analyzing puzzle...";
            NotifyVisualChanged();
            try
            {
                await Task.Delay(20);
                return AnalyzeCurrentPuzzle();
            }
            finally
            {
                IsBusy = false;
                BusyMessage = null;
                NotifyVisualChanged();
            }
        }

        public bool ToggleEdge(int edgeIndex, bool isAlternative)
        {
            if (CurrentMesh == null || IsGenerating || IsBusy || edgeIndex < 0 || edgeIndex >= CurrentMesh.Edges.Count)
                return false;

            PendingGroupEdgeIndex = null;
            bool dismissedVisual = DismissVisualSolveInternal();

            if (MarkedEdges.Contains(edgeIndex))
            {
                if (dismissedVisual) NotifyVisualChanged();
                return false;
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
                return true;
            }

            if (dismissedVisual) NotifyVisualChanged();
            if (DisallowFalseMove)
                FlashInvalidEdge(edgeIndex);
            return false;
        }

        public bool SetEdgeDirect(int edgeIndex, EdgeState targetState)
        {
            if (CurrentMesh == null || IsGenerating || IsBusy || edgeIndex < 0 || edgeIndex >= CurrentMesh.Edges.Count)
                return false;

            bool dismissedVisual = DismissVisualSolveInternal();

            if (MarkedEdges.Contains(edgeIndex) || CurrentMesh.Edges[edgeIndex].State == targetState)
            {
                if (dismissedVisual) NotifyVisualChanged();
                return false;
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
                return true;
            }

            if (dismissedVisual) NotifyVisualChanged();
            if (DisallowFalseMove)
                FlashInvalidEdge(edgeIndex);
            return false;
        }

        public bool ToggleCellColor(int cellIndex, bool isAlternative)
        {
            if (CurrentMesh == null || IsGenerating || IsBusy || cellIndex < 0 || cellIndex >= CurrentMesh.Cells.Count)
                return false;

            bool dismissedVisual = DismissVisualSolveInternal();
            var action = new PuzzleCellColorAction(CurrentMesh, cellIndex, isAlternative);
            if (UndoTree.Do(action))
            {
                NotifyChanged();
                return true;
            }
            if (dismissedVisual) NotifyVisualChanged();
            return false;
        }

        public void Undo()
        {
            if (IsGenerating || IsBusy) return;
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
                NotifyVisualChanged();
            }
        }

        public bool UndoAndForget()
        {
            if (IsGenerating || IsBusy) return false;
            bool dismissedVisual = DismissVisualSolveInternal();
            if (UndoTree.CanUndo && UndoTree.UndoAndForget())
            {
                WasAutoSolved = false;
                CheckSolution();
                NotifyChanged();
                return true;
            }
            if (dismissedVisual) NotifyVisualChanged();
            return false;
        }

        public void Redo()
        {
            if (IsGenerating || IsBusy) return;
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
                NotifyVisualChanged();
            }
        }

        public void Fix()
        {
            if (CurrentMesh == null || IsGenerating || IsBusy) return;
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
            if (IsGenerating || IsBusy) return;
            DismissVisualSolveInternal();
            MarkedEdges.Clear();
            UndoTree.ClearMark();
            preVictoryMarkedEdges = null;
            preVictoryUndoMark = null;
            preVictoryUndoTree = null;
            CheckSolution();
            NotifyChanged();
        }

        public void RevertToFix()
        {
            if (IsGenerating || IsBusy) return;
            bool dismissedVisual = DismissVisualSolveInternal();
            if (!UndoTree.CanRevert)
            {
                if (dismissedVisual) NotifyVisualChanged();
                return;
            }
            UndoTree.RevertToMark();
            WasAutoSolved = UndoTree.CurrentAction is PuzzleBatchSolveAction;
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
                if (StartTime.HasValue)
                {
                    solvedElapsed = DateTime.UtcNow - StartTime.Value;
                }
                preVictoryMarkedEdges = MarkedEdges.ToArray();
                preVictoryUndoMark = UndoTree.GetMarkedDirect();
                preVictoryUndoTree = UndoTree;
                IsSolved = true;
                Fix(); // Auto-fix solution upon victory
            }
            else if (!currentlySolved && IsSolved)
            {
                IsSolved = false;
                if (solvedElapsed.HasValue)
                {
                    StartTime = DateTime.UtcNow - solvedElapsed.Value;
                    solvedElapsed = null;
                }
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

            string puzzleSerialized = (ActiveWidth > 0 && ActiveHeight > 0)
                ? PuzzleCodec.Encode(CurrentMesh, ActiveWidth, ActiveHeight)
                : (PuzzleFormatHelper.TryGetSquareDimensions(CurrentMesh, out int sw, out int sh)
                    ? PuzzleCodec.Encode(CurrentMesh, sw, sh)
                    : CurrentMesh.SaveToString());

            string activeSize = (ActiveWidth > 0 && ActiveHeight > 0)
                ? $"{ActiveWidth}x{ActiveHeight}"
                : SizeText;

            return new SavedGameState
            {
                Puzzle = puzzleSerialized,
                EdgeStates = edgeStates,
                CellColors = cellColors,
                MarkedEdges = MarkedEdges.ToArray(),
                ElapsedSeconds = ElapsedTime.TotalSeconds,
                IsSolved = IsSolved,
                WasAutoSolved = WasAutoSolved,
                DisallowFalseMove = DisallowFalseMove,
                ShowCellColors = ShowCellColors,
                Type = ActiveType.ToString(),
                Size = activeSize,
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
            width = 0;
            height = 0;
            try
            {
                var lines = data.Split(LineSeparators, StringSplitOptions.RemoveEmptyEntries);
                mesh = new Mesh(0, 0, MeshType.Square);
                if (!mesh.LoadFromText(lines))
                {
                    mesh = null;
                    return false;
                }
                if (PuzzleFormatHelper.TryGetSquareDimensions(mesh, out int sw, out int sh))
                {
                    width = sw;
                    height = sh;
                }
                return true;
            }
            catch
            {
                mesh = null;
                return false;
            }
        }

        public bool RestoreState(SavedGameState state)
        {
            if (string.IsNullOrEmpty(state.Puzzle)) return false;

            if (!TryDecodeOrLoadMesh(state.Puzzle, out var mesh, out bool isDecoded, out int w, out int h) || mesh == null)
                return false;

            CurrentMesh = mesh;
            CurrentType = mesh.MeshType;
            if (w > 0 && h > 0)
            {
                ActiveWidth = w;
                ActiveHeight = h;
                SizeText = $"{w}x{h}";
            }
            else if (!string.IsNullOrEmpty(state.Size) && PuzzleHelper.ParseSize(state.Size, mesh.MeshType, out int pw, out int ph))
            {
                ActiveWidth = pw;
                ActiveHeight = ph;
                SizeText = state.Size;
            }
            else
            {
                ActiveWidth = 0;
                ActiveHeight = 0;
            }

            Difficulty = Math.Clamp(state.Difficulty, 0, 5);
            TotalIntersections = mesh.Intersections.Count;
            PrunedProgress = 0;
            IsGenerating = false;
            UndoTree = new UndoTree();
            MarkedEdges.Clear();

            // Replay edge states and cell colors through Mesh actions so edgeSets and cellColorSets stay consistent
            if (state.EdgeStates != null && state.EdgeStates.Length == mesh.Edges.Count)
            {
                mesh.Clear();
                for (int i = 0; i < state.EdgeStates.Length; i++)
                {
                    EdgeState es = (EdgeState)state.EdgeStates[i];
                    if (Enum.IsDefined(es) && es != EdgeState.Empty)
                    {
                        new SetAction(mesh, i, es).Perform();
                    }
                }
            }
            mesh.RecalculateCounts();

            if (state.CellColors != null && state.CellColors.Length == mesh.Cells.Count)
            {
                for (int i = 0; i < state.CellColors.Length; i++)
                {
                    int color = state.CellColors[i];
                    if (color is 1 or -1)
                    {
                        new CellColorJoinAction(mesh, i, -1, color == 1).Perform();
                    }
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

            double elapsedSec = double.IsFinite(state.ElapsedSeconds) && state.ElapsedSeconds > 0
                ? state.ElapsedSeconds
                : 0;
            StartTime = elapsedSec > 0
                ? DateTime.UtcNow.AddSeconds(-elapsedSec)
                : DateTime.UtcNow;

            IsSolved = state.IsSolved;
            solvedElapsed = state.IsSolved ? TimeSpan.FromSeconds(elapsedSec) : null;
            WasAutoSolved = state.IsSolved && state.WasAutoSolved;
            ClearTransientBoardState();
            CheckSolution();
            NotifyVisualChanged();
            return true;
        }

        public bool LoadPuzzleText(string puzzleText)
        {
            if (string.IsNullOrEmpty(puzzleText)) return false;

            if (!TryDecodeOrLoadMesh(puzzleText, out var mesh, out bool isDecoded, out int w, out int h) || mesh == null)
                return false;

            ActiveWidth = w;
            ActiveHeight = h;
            if (w > 0 && h > 0)
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
            solvedElapsed = null;
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
            NotifyVisualChanged();

            _ = ClearAfterDelayAsync(600, () =>
            {
                if (InvalidEdgeIndex == edgeIndex)
                {
                    InvalidEdgeIndex = null;
                    NotifyVisualChanged();
                }
            }, token);
        }

        public void ClearHint()
        {
            ClearHintInternal();
            NotifyVisualChanged();
        }

        public void ResetCurrentPuzzle()
        {
            if (CurrentMesh == null || IsGenerating || IsBusy) return;

            CurrentMesh.Clear();
            MarkedEdges.Clear();
            UndoTree = new UndoTree();
            ApplyAutoStartIfEnabled(CurrentMesh);
            StartTime = DateTime.UtcNow;
            solvedElapsed = null;
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
                NotifyVisualChanged();
        }

        public bool StartVisualSolve()
        {
            if (CurrentMesh == null || IsGenerating) return false;

            DismissVisualSolveInternal();
            ClearHintInternal();

            Mesh solver = PuzzleHelper.CreateFullSolverMesh(CurrentMesh);
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
                NotifyVisualChanged();
                return true;
            }

            DismissVisualSolveInternal();
            NotifyVisualChanged();
            return false;
        }

        public async Task<bool> StartVisualSolveAsync()
        {
            if (CurrentMesh == null || IsGenerating || IsBusy) return false;
            IsBusy = true;
            BusyMessage = "Preparing visual solver...";
            NotifyVisualChanged();
            try
            {
                await Task.Delay(20);
                return StartVisualSolve();
            }
            finally
            {
                IsBusy = false;
                BusyMessage = null;
                NotifyVisualChanged();
            }
        }

        public void SetVisualSolveStep(int index)
        {
            if (!IsVisualSolving) return;
            int clamped = Math.Clamp(index, 0, visualSolveSteps.Count - 1);
            if (clamped == VisualSolveStepIndex) return;
            VisualSolveStepIndex = clamped;
            NotifyVisualChanged();
        }

        public bool StepVisualSolveForward()
        {
            if (!IsVisualSolving || VisualSolveStepIndex >= visualSolveSteps.Count - 1)
                return false;
            VisualSolveStepIndex++;
            NotifyVisualChanged();
            return true;
        }

        public bool StepVisualSolveBackward()
        {
            if (!IsVisualSolving || VisualSolveStepIndex <= 0)
                return false;
            VisualSolveStepIndex--;
            NotifyVisualChanged();
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
            Mesh solver = PuzzleHelper.CreateFullSolverMesh(CurrentMesh);
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

        public async Task<bool> SolvePuzzleAsync()
        {
            if (CurrentMesh == null || IsGenerating || IsBusy) return false;
            IsBusy = true;
            BusyMessage = "Solving puzzle...";
            NotifyVisualChanged();
            try
            {
                await Task.Delay(20);
                return SolvePuzzle();
            }
            finally
            {
                IsBusy = false;
                BusyMessage = null;
                NotifyVisualChanged();
            }
        }

        public bool RequestHint(bool applyMove)
        {
            if (CurrentMesh == null || IsGenerating) return false;
            DismissVisualSolveInternal();

            var (edgeIdx, targetState, message) = ComputeHint();
            HintMessage = message;
            bool boardMutated = false;

            if (edgeIdx.HasValue)
            {
                HintEdgeIndex = edgeIdx.Value;
                HintTargetState = targetState;
                HintIsApplied = false;
                if (applyMove)
                {
                    var action = new PuzzleSetEdgeStateAction(CurrentMesh, edgeIdx.Value, targetState, DisallowFalseMove);
                    if (UndoTree.Do(action))
                    {
                        HintIsApplied = true;
                        boardMutated = true;
                        CheckSolution();
                    }
                }
                else
                {
                    // Don't give away whether it's a line or a cross.
                    HintMessage = "The highlighted edge can be determined from what's around it.";
                }

                var flashToken = ResetTokenSource(ref hintFlashCts);
                int delayMs = HintIsApplied ? 1500 : 3500;
                int capturedEdge = edgeIdx.Value;
                _ = ClearAfterDelayAsync(delayMs, () =>
                {
                    if (HintEdgeIndex == capturedEdge)
                    {
                        HintEdgeIndex = null;
                        NotifyVisualChanged();
                    }
                }, flashToken);
            }

            var msgToken = ResetTokenSource(ref hintMessageCts);
            _ = ClearAfterDelayAsync(5000, () =>
            {
                HintMessage = null;
                NotifyVisualChanged();
            }, msgToken);

            if (boardMutated)
                NotifyChanged();
            else
                NotifyVisualChanged();
            return edgeIdx.HasValue;
        }

        public async Task<bool> RequestHintAsync(bool applyMove)
        {
            if (CurrentMesh == null || IsGenerating || IsBusy) return false;
            IsBusy = true;
            BusyMessage = "Computing hint...";
            NotifyVisualChanged();
            try
            {
                await Task.Delay(20);
                return RequestHint(applyMove);
            }
            finally
            {
                IsBusy = false;
                BusyMessage = null;
                NotifyVisualChanged();
            }
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
            Mesh fullSolverMesh = PuzzleHelper.CreateFullSolverMesh(CurrentMesh);
            SolveState solveRes = fullSolverMesh.TrySolve();
            if ((solveRes == SolveState.Solved || solveRes == SolveState.MultipleSolutions) && fullSolverMesh.SolutionFound != null)
            {
                var sol = fullSolverMesh.SolutionFound;
                for (int i = 0; i < CurrentMesh.Edges.Count; i++)
                {
                    if (CurrentMesh.Edges[i].State == EdgeState.Empty && sol.Edges[i].State != EdgeState.Empty)
                    {
                        string stateStr = sol.Edges[i].State == EdgeState.Filled ? "line" : "cross (✕)";
                        string prefix = solveRes == SolveState.MultipleSolutions
                            ? "Suggested move (multiple solutions exist): This edge can be a "
                            : "Advanced deduction: This edge must be a ";
                        return (i, sol.Edges[i].State, $"{prefix}{stateStr}.");
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
            if (generatingMesh != null)
            {
                generatingMesh.AbortPrune = true;
            }
            CancelAndDispose(ref generationCts);
            CancelAndDispose(ref invalidFlashCts);
            CancelAndDispose(ref hintFlashCts);
            CancelAndDispose(ref hintMessageCts);
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
