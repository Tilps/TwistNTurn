using System.Text.Json;
using TwistNTurn.Web.Services;

[assembly: Parallelize(Scope = ExecutionScope.MethodLevel)]

namespace TwistNTurn.Tests
{
    [TestClass]
    public class QualityRegressionTests
    {
        private static Mesh CreateSquare2x2CornerBlackPuzzle()
        {
            var mesh = new Mesh(2, 2, MeshType.Square);
            mesh.Intersections[0].Type = IntersType.Turn;
            mesh.Intersections[2].Type = IntersType.Turn;
            mesh.Intersections[6].Type = IntersType.Turn;
            mesh.Intersections[8].Type = IntersType.Turn;
            return mesh;
        }

        [TestMethod]
        public void PerformUnsetZero_SplitsEdgeSetsCorrectly_AndGuardsEmptyEdge()
        {
            var mesh = new Mesh(3, 1, MeshType.Square);
            int e0 = mesh.GetEdgeJoining(0, 2);
            int e1 = mesh.GetEdgeJoining(2, 4);
            int e2 = mesh.GetEdgeJoining(4, 6);

            Assert.IsTrue(new SetAction(mesh, e0, EdgeState.Filled).Perform());
            Assert.IsTrue(new SetAction(mesh, e1, EdgeState.Filled).Perform());
            Assert.IsTrue(new SetAction(mesh, e2, EdgeState.Filled).Perform());

            Assert.AreNotEqual(0, mesh.Edges[e0].EdgeSet);
            Assert.AreEqual(mesh.Edges[e0].EdgeSet, mesh.Edges[e1].EdgeSet);
            Assert.AreEqual(mesh.Edges[e0].EdgeSet, mesh.Edges[e2].EdgeSet);

            // Break middle edge e1: e0 and e2 must belong to distinct non-zero EdgeSets
            var unsetMiddle = new UnsetAction(mesh, e1);
            Assert.IsTrue(unsetMiddle.Perform());
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[e1].State);
            Assert.AreEqual(0, mesh.Edges[e1].EdgeSet);
            Assert.AreNotEqual(0, mesh.Edges[e0].EdgeSet);
            Assert.AreNotEqual(0, mesh.Edges[e2].EdgeSet);
            Assert.AreNotEqual(mesh.Edges[e0].EdgeSet, mesh.Edges[e2].EdgeSet,
                "Disconnected filled chains must have distinct EdgeSet IDs after unset");

            // Unperforming the unset must rejoin them into the original EdgeSet
            unsetMiddle.Unperform();
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[e1].State);
            Assert.AreEqual(mesh.Edges[e0].EdgeSet, mesh.Edges[e2].EdgeSet);

            // Unsetting an already-empty edge must fail cleanly without corrupting ExcludedCount
            int emptyEdge = mesh.GetEdgeJoining(0, 1);
            var unsetEmpty = new UnsetAction(mesh, emptyEdge);
            Assert.IsFalse(unsetEmpty.Perform());
            Assert.IsFalse(unsetEmpty.Successful);
            Assert.AreEqual(0, mesh.Cells[0].ExcludedCount);
            Assert.AreEqual(0, mesh.Intersections[0].ExcludedCount);
        }

        [TestMethod]
        public void MeshLoadFromText_PreservesNonSquareMeshType_AndRejectsMalformedInputs()
        {
            var hex2 = new Mesh(5, 5, MeshType.Hexagonal2);
            hex2.Intersections[0].Type = IntersType.Twist;
            hex2.Intersections[1].Type = IntersType.Turn;
            string serialized = hex2.SaveToString();

            Assert.IsTrue(PuzzleFormatHelper.TryParsePuzzle(serialized, out Mesh? loaded, out string fmt, out _, out _));
            Assert.AreEqual("Twist File", fmt);
            Assert.IsNotNull(loaded);
            Assert.AreEqual(MeshType.Hexagonal2, loaded.MeshType, "LoadFromText must preserve non-square MeshType from line 0");
            Assert.HasCount(hex2.Intersections.Count, loaded.Intersections);
            Assert.AreEqual(IntersType.Twist, loaded.Intersections[0].Type);
            Assert.AreEqual(IntersType.Turn, loaded.Intersections[1].Type);

            string[] badPayloads =
            [
                "Square\nIntersections\n2\n0 0 Unknown", // truncated
                "Square\nIntersections\n1\nNaN 0 Unknown\nEdges\n0", // NaN coordinate
                "Square\nIntersections\n2\n0 0 Unknown\n1 0 Unknown\nEdges\n1\n0 99 Filled", // out-of-range vertex
                "Square\nIntersections\n2\n0 0 Unknown\n1 0 Unknown\nEdges\n1\n0 0 Filled", // self-loop
                "Square\nIntersections\n2\n0 0 Unknown\n1 0 Unknown\nEdges\n2\n0 1 Filled\n1 0 Excluded", // duplicate edge
                "NotARealMeshType\nIntersections\n0\nEdges\n0" // invalid MeshType
            ];

            foreach (string payload in badPayloads)
            {
                var target = new Mesh(0, 0, MeshType.Square);
                Assert.IsFalse(target.LoadFromText(payload.AsSpan()), $"Expected malformed payload to fail cleanly: {payload}");
            }
        }

        [TestMethod]
        public void PuzzleCodec_RejectsArbitraryUrls_OversizedGrids_AndInvalidClues()
        {
            Assert.IsFalse(PuzzleCodec.TryDecode("https://example.com/foo", out _, out _, out _, out _));
            Assert.IsFalse(PuzzleCodec.TryDecode("unknown/3x3/", out _, out _, out _, out _));
            Assert.IsFalse(PuzzleCodec.TryDecode("sq/500x500/", out _, out _, out _, out _));
            Assert.IsFalse(PuzzleCodec.TryDecode("sq//", out _, out _, out _, out _));
            Assert.IsFalse(PuzzleCodec.TryDecode("sq/2x2/2!0", out _, out _, out _, out _));
            Assert.IsFalse(PuzzleCodec.TryDecode("500x500:0", out _, out _, out _, out _));
            Assert.IsFalse(PuzzleHelper.ParseSize("1000x1000", MeshType.Square, out _, out _));
            Assert.AreEqual(MeshType.Square, PuzzleHelper.MeshTypeFromString("99999"));
            Assert.AreEqual(MeshType.Hexagonal4, PuzzleHelper.MeshTypeFromString("hexagon 4"));
        }

        [TestMethod]
        public void PuzzleCompression_RejectsOversizedDeflatePayload()
        {
            string bigText = new string('A', 32 * 1024);
            string compressed = PuzzleCompression.CompressToUrlSafe(bigText);

            Assert.IsFalse(PuzzleCompression.TryDecompressFromUrlSafe(compressed, out _, maxBytes: 1024));
            Assert.IsTrue(PuzzleCompression.TryDecompressFromUrlSafe(compressed, out string ok, maxBytes: 64 * 1024));
            Assert.AreEqual(bigText, ok);
        }

        [TestMethod]
        public void PuzzleFormatHelper_StrictAscii_AndPreservesBlankVerticalLines()
        {
            // Bogus text starting with '+' must be rejected
            Assert.IsFalse(PuzzleFormatHelper.TryParsePuzzle("+abc+\nhello\n+xyz+", out _, out _, out _, out _));

            // Unsolved ASCII grid whose empty vertical-edge lines had trailing whitespace trimmed must still parse
            var mesh = CreateSquare2x2CornerBlackPuzzle();
            string? ascii = PuzzleFormatHelper.ExportSquareAscii(mesh);
            Assert.IsNotNull(ascii);
            string trimmedLines = string.Join("\n", ascii.Replace("\r\n", "\n").Replace("\r", "\n").Split('\n').Select(l => l.TrimEnd()));
            Assert.IsTrue(PuzzleFormatHelper.TryParsePuzzle(trimmedLines, out Mesh? parsed, out string fmt, out int w, out int h));
            Assert.AreEqual("ASCII Grid", fmt);
            Assert.AreEqual(2, w);
            Assert.AreEqual(2, h);
            Assert.IsNotNull(parsed);
            Assert.AreEqual(IntersType.Turn, parsed.Intersections[0].Type);
            Assert.AreEqual(IntersType.Turn, parsed.Intersections[8].Type);
        }

        [TestMethod]
        public void UndoTree_RevertToMark_IgnoresNullOrDetachedMark_AndUndoAndForgetUpdatesMark()
        {
            var mesh = new Mesh(2, 2, MeshType.Square);
            var tree = new UndoTree();

            Assert.IsTrue(tree.Do(new PuzzleSetEdgeStateAction(mesh, 0, EdgeState.Filled)));
            Assert.IsTrue(tree.Do(new PuzzleSetEdgeStateAction(mesh, 1, EdgeState.Filled)));

            // RevertToMark with no mark must be a no-op
            tree.RevertToMark();
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[1].State);

            // Mark current node, then UndoAndForget it -> mark moves cleanly to parent
            tree.Mark();
            Assert.IsTrue(tree.UndoAndForget());
            Assert.IsTrue(tree.Do(new PuzzleSetEdgeStateAction(mesh, 2, EdgeState.Filled)));
            Assert.IsTrue(tree.CanRevert);
            tree.RevertToMark();
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[2].State);
        }

        [TestMethod]
        public void GameService_PendingGridAndSizeEdits_DoNotCorruptActiveExportOrShare()
        {
            using var game = new GameService();
            var sq2 = CreateSquare2x2CornerBlackPuzzle();
            string expectedCode = PuzzleCodec.Encode(sq2, 2, 2);
            Assert.IsTrue(game.LoadPuzzleText(expectedCode));
            Assert.AreEqual(2, game.ActiveWidth);
            Assert.AreEqual(2, game.ActiveHeight);

            // User edits pending form inputs without clicking "New Puzzle"
            game.CurrentType = MeshType.Hexagonal2;
            game.SizeText = "7x9";

            Assert.AreEqual(expectedCode, game.GetShareCode());

            var exported = game.ExportState();
            Assert.IsNotNull(exported);
            Assert.AreEqual(expectedCode, exported.Puzzle);
            Assert.AreEqual("2x2", exported.Size);
            Assert.AreEqual("Square", exported.Type);
        }

        [TestMethod]
        public void GameService_RestoreState_ReplaysCellColorsSafely_PreservesPreferences_AndFreezesElapsed()
        {
            using var game = new GameService();
            game.DisallowFalseMove = true;
            game.ShowCellColors = true;

            string code = PuzzleCodec.Encode(CreateSquare2x2CornerBlackPuzzle(), 2, 2);
            var state = new SavedGameState
            {
                Puzzle = code,
                EdgeStates = [99, -5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                CellColors = [1, 42, -1, 0],
                MarkedEdges = [-1, 0, 999],
                ElapsedSeconds = double.NaN,
                IsSolved = false,
                DisallowFalseMove = false,
                ShowCellColors = false
            };

            Assert.IsTrue(game.RestoreState(state));
            Assert.IsTrue(game.DisallowFalseMove, "RestoreState must not overwrite UserPreferences.DisallowFalseMove");
            Assert.IsTrue(game.ShowCellColors, "RestoreState must not overwrite UserPreferences.ShowCellColors");
            Assert.AreEqual(EdgeState.Empty, game.CurrentMesh!.Edges[0].State);
            Assert.AreEqual(EdgeState.Empty, game.CurrentMesh.Edges[1].State);
            Assert.AreEqual(EdgeState.Filled, game.CurrentMesh.Edges[2].State);
            Assert.AreEqual(1, game.CurrentMesh.Cells[0].Color);
            Assert.AreEqual(0, game.CurrentMesh.Cells[1].Color);
            Assert.AreEqual(-1, game.CurrentMesh.Cells[2].Color);

            // Clicking a restored shaded cell must NOT crash with ArgumentOutOfRangeException
            Assert.IsTrue(game.ToggleCellColor(0, isAlternative: false));
            Assert.AreEqual(-1, game.CurrentMesh.Cells[0].Color);
            Assert.IsTrue(game.ToggleCellColor(2, isAlternative: false));
            Assert.AreEqual(0, game.CurrentMesh.Cells[2].Color);

            // Solving freezes ElapsedTime
            Assert.IsTrue(game.SolvePuzzle());
            Assert.IsTrue(game.IsSolved);
            TimeSpan frozen1 = game.ElapsedTime;
            Thread.Sleep(25);
            TimeSpan frozen2 = game.ElapsedTime;
            Assert.AreEqual(frozen1, frozen2, "ElapsedTime must stay frozen once solved");
        }

        [TestMethod]
        public void AotJsonContext_RoundTripsSavedGameState()
        {
            var saved = new SavedGameState
            {
                Puzzle = PuzzleCodec.Encode(CreateSquare2x2CornerBlackPuzzle(), 2, 2),
                EdgeStates = [1, 2, 0],
                CellColors = [1, -1, 0, 0],
                MarkedEdges = [0],
                ElapsedSeconds = 123.5,
                IsSolved = true,
                WasAutoSolved = false,
                Type = "Square",
                Size = "2x2",
                Difficulty = 2
            };
            string json = JsonSerializer.Serialize(saved, TwistJsonContext.Default.SavedGameState);
            var back = JsonSerializer.Deserialize(json, TwistJsonContext.Default.SavedGameState);
            Assert.IsNotNull(back);
            Assert.AreEqual(saved.Puzzle, back.Puzzle);
            Assert.AreEqual(saved.ElapsedSeconds, back.ElapsedSeconds);
            Assert.IsTrue(back.IsSolved);
        }

        [TestMethod]
        public void MeshCopyConstructor_CopiesAllOptions_SeededGenerationIsDeterministic_AndPercentSolvedIsNotInverted()
        {
            var src = new Mesh(2, 2, MeshType.Square)
            {
                UseMerging = true,
                UseCellColoring = true,
                UseCellColoringTrials = false,
                UseDerivedColoring = false,
                SuperSlowMo = true,
                GenerateBoringFraction = 0.05
            };

            var clone = new Mesh(src);
            Assert.IsTrue(clone.UseMerging);
            Assert.IsTrue(clone.UseCellColoring);
            Assert.IsFalse(clone.UseCellColoringTrials);
            Assert.IsFalse(clone.UseDerivedColoring);
            Assert.IsTrue(clone.SuperSlowMo);
            Assert.AreEqual(0.05, clone.GenerateBoringFraction, 1e-6);

            // Deterministic seeded generation produces identical puzzles
            var m1 = PuzzleHelper.MakeMesh(3, 3, MeshType.Square, 1);
            var m2 = PuzzleHelper.MakeMesh(3, 3, MeshType.Square, 1);
            m1.Generate(new Random(12345));
            m2.Generate(new Random(12345));
            Assert.AreEqual(PuzzleCodec.Encode(m1, 3, 3), PuzzleCodec.Encode(m2, 3, 3));

            // PercentSolved on iterative solve must be 1.0 when solved
            var solvedMesh = CreateSquare2x2CornerBlackPuzzle();
            solvedMesh.SetRatingCodeOptions("S");
            Assert.AreEqual(SolveState.Solved, solvedMesh.TrySolve());
            Assert.AreEqual(1.0, solvedMesh.PercentSolved, 1e-6);
        }

        [TestMethod]
        public void ToggleEdge_ReturnsFalseWhenRejectedByDisallowFalseMove()
        {
            using var game = new GameService();
            var mesh = new Mesh(2, 2, MeshType.Square);
            int center = 4;
            var centerEdges = mesh.Intersections[center].Edges;
            new SetAction(mesh, centerEdges[0], EdgeState.Filled).Perform();
            new SetAction(mesh, centerEdges[1], EdgeState.Filled).Perform();
            game.LoadPuzzleText(PuzzleFormatHelper.ExportTwistText(mesh));
            game.DisallowFalseMove = true;

            Assert.IsFalse(game.ToggleEdge(centerEdges[2], isAlternative: false));
            Assert.IsFalse(game.CanUndo);
            Assert.AreEqual(centerEdges[2], game.InvalidEdgeIndex);
        }
    }
}
