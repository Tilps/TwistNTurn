using Microsoft.VisualStudio.TestTools.UnitTesting;
using TwistNTurn;
using TwistNTurn.Web.Services;

namespace TwistNTurn.Tests
{
    [TestClass]
    public class Phase1FeaturesTests
    {
        private static Mesh CreateSquare2x2CornerBlackPuzzle()
        {
            var mesh = new Mesh(2, 2, MeshType.Square);
            // 2x2 cell grid -> 3x3 intersections: corners are 0=(0,0), 2=(0,2), 6=(2,0), 8=(2,2)
            mesh.Intersections[0].Type = IntersType.Turn;
            mesh.Intersections[2].Type = IntersType.Turn;
            mesh.Intersections[6].Type = IntersType.Turn;
            mesh.Intersections[8].Type = IntersType.Turn;
            return mesh;
        }

        [TestMethod]
        public void ComputeHint_DetectsMasyuContradictions()
        {
            var game = new GameService();

            // 1. Intersection degree > 2
            var degMesh = new Mesh(2, 2, MeshType.Square);
            int center = 4; // (1,1) in 3x3 intersections
            var centerEdges = degMesh.Intersections[center].Edges;
            new SetAction(degMesh, centerEdges[0], EdgeState.Filled).Perform();
            new SetAction(degMesh, centerEdges[1], EdgeState.Filled).Perform();
            new SetAction(degMesh, centerEdges[2], EdgeState.Filled).Perform();
            game.LoadPuzzleText(PuzzleFormatHelper.ExportTwistText(degMesh));

            var degHint = game.ComputeHint();
            Assert.IsNull(degHint.EdgeIndex);
            Assert.Contains("more than 2", degHint.Message);

            // 2. White pearl (Twist) with a turn
            var twistTurnMesh = new Mesh(2, 2, MeshType.Square);
            twistTurnMesh.Intersections[4].Type = IntersType.Twist; // White pearl at (1,1)
            int topEdge = twistTurnMesh.GetEdgeJoining(4, 3);      // (1,1)-(1,0)
            int rightEdge = twistTurnMesh.GetEdgeJoining(4, 7);    // (1,1)-(2,1)
            new SetAction(twistTurnMesh, topEdge, EdgeState.Filled).Perform();
            new SetAction(twistTurnMesh, rightEdge, EdgeState.Filled).Perform();
            game.LoadPuzzleText(PuzzleFormatHelper.ExportTwistText(twistTurnMesh));

            var twistHint = game.ComputeHint();
            Assert.IsNull(twistHint.EdgeIndex);
            Assert.Contains("white pearl", twistHint.Message);

            // 3. Black pearl (Turn) going straight through
            var turnStraightMesh = new Mesh(2, 2, MeshType.Square);
            turnStraightMesh.Intersections[4].Type = IntersType.Turn; // Black pearl at (1,1)
            int topE = turnStraightMesh.GetEdgeJoining(4, 3);         // (1,1)-(1,0)
            int bottomE = turnStraightMesh.GetEdgeJoining(4, 5);      // (1,1)-(1,2)
            new SetAction(turnStraightMesh, topE, EdgeState.Filled).Perform();
            new SetAction(turnStraightMesh, bottomE, EdgeState.Filled).Perform();
            game.LoadPuzzleText(PuzzleFormatHelper.ExportTwistText(turnStraightMesh));

            var turnHint = game.ComputeHint();
            Assert.IsNull(turnHint.EdgeIndex);
            Assert.Contains("black pearl", turnHint.Message);
        }

        [TestMethod]
        public void RequestHint_SeparatesHighlightFromApply_AndSupportsUndo()
        {
            var game = new GameService();
            var mesh = CreateSquare2x2CornerBlackPuzzle();
            game.LoadPuzzleText(PuzzleCodec.Encode(mesh, 2, 2));

            // Plain Hint highlights without modifying board
            Assert.IsTrue(game.RequestHint(applyMove: false));
            Assert.IsNotNull(game.HintEdgeIndex);
            Assert.IsFalse(game.HintIsApplied);
            Assert.IsFalse(game.CanUndo);

            // Do Hint applies a logical edge and supports Undo
            Assert.IsTrue(game.RequestHint(applyMove: true));
            Assert.IsTrue(game.HintIsApplied);
            Assert.IsNotNull(game.HintEdgeIndex);
            int appliedEdge = game.HintEdgeIndex.Value;
            Assert.AreNotEqual(EdgeState.Empty, game.CurrentMesh!.Edges[appliedEdge].State);
            Assert.IsTrue(game.CanUndo);

            game.Undo();
            Assert.AreEqual(EdgeState.Empty, game.CurrentMesh.Edges[appliedEdge].State);
        }

        [TestMethod]
        public void SolvePuzzle_AndReset_WorkWithUndoAndVictoryLockRestore()
        {
            var game = new GameService();
            var mesh = CreateSquare2x2CornerBlackPuzzle();
            game.LoadPuzzleText(PuzzleCodec.Encode(mesh, 2, 2));

            // Lock one unrelated move before solving so we can verify pre-victory lock restore on Undo
            game.Fix();
            int preFixCount = game.MarkedEdges.Count;

            Assert.IsTrue(game.SolvePuzzle());
            Assert.IsTrue(game.IsSolved);
            Assert.IsTrue(game.WasAutoSolved);
            Assert.IsTrue(game.MarkedEdges.Count > preFixCount);

            // Undo reverses the solve AND restores pre-victory lock state
            game.Undo();
            Assert.IsFalse(game.IsSolved);
            Assert.IsFalse(game.WasAutoSolved);
            Assert.HasCount(preFixCount, game.MarkedEdges);

            // Set an edge then Reset
            game.SetEdgeDirect(0, EdgeState.Filled);
            Assert.IsTrue(game.CanUndo);
            game.ResetCurrentPuzzle();
            Assert.IsFalse(game.CanUndo);
            Assert.AreEqual(EdgeState.Empty, game.CurrentMesh!.Edges[0].State);
        }

        [TestMethod]
        public void DisallowFalseMove_FlashesInvalidEdgeOnRejection()
        {
            var game = new GameService();
            var mesh = new Mesh(2, 2, MeshType.Square);
            int center = 4; // (1,1)
            var edges = mesh.Intersections[center].Edges;
            new SetAction(mesh, edges[0], EdgeState.Filled).Perform();
            new SetAction(mesh, edges[1], EdgeState.Filled).Perform();
            game.LoadPuzzleText(PuzzleFormatHelper.ExportTwistText(mesh));
            game.DisallowFalseMove = true;

            game.SetEdgeDirect(edges[2], EdgeState.Filled);
            Assert.AreEqual(EdgeState.Empty, game.CurrentMesh!.Edges[edges[2]].State);
            Assert.AreEqual(edges[2], game.InvalidEdgeIndex);
        }
    }
}
