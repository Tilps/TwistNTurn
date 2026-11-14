using Microsoft.VisualStudio.TestTools.UnitTesting;
using TwistNTurn;
using TwistNTurn.Web.Services;

namespace TwistNTurn.Tests
{
    [TestClass]
    public class Phase4And5Tests
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
        public void PuzzleFormatHelper_RoundTripsTwistText_Ascii_PearlDescriptor_And_ShareLink()
        {
            var mesh = new Mesh(3, 2, MeshType.Square);
            mesh.Intersections[0].Type = IntersType.Twist; // White pearl at (0,0)
            mesh.Intersections[5].Type = IntersType.Turn;  // Black pearl at (1,2)
            new SetAction(mesh, 0, EdgeState.Filled).Perform();
            new SetAction(mesh, 1, EdgeState.Excluded).Perform();

            Assert.IsTrue(PuzzleFormatHelper.CanExportSquareFormats(mesh));
            Assert.IsTrue(PuzzleFormatHelper.TryGetSquareDimensions(mesh, out int w, out int h));
            Assert.AreEqual(3, w);
            Assert.AreEqual(2, h);

            var game = new GameService();

            // 1. Native .twist format
            string twistText = PuzzleFormatHelper.ExportTwistText(mesh);
            Assert.IsTrue(game.ImportPuzzle(twistText, out string twistFmt));
            Assert.AreEqual("Twist File", twistFmt);
            Assert.AreEqual(IntersType.Twist, game.CurrentMesh!.Intersections[0].Type);
            Assert.AreEqual(IntersType.Turn, game.CurrentMesh.Intersections[5].Type);
            Assert.AreEqual(EdgeState.Filled, game.CurrentMesh.Edges[0].State);
            Assert.AreEqual(EdgeState.Excluded, game.CurrentMesh.Edges[1].State);

            // 2. Square Masyu ASCII format
            string? ascii = PuzzleFormatHelper.ExportSquareAscii(mesh);
            Assert.IsNotNull(ascii);
            Assert.IsTrue(game.ImportPuzzle(ascii, out string asciiFmt));
            Assert.AreEqual("ASCII Grid", asciiFmt);
            Assert.AreEqual(IntersType.Twist, game.CurrentMesh!.Intersections[0].Type);
            Assert.AreEqual(IntersType.Turn, game.CurrentMesh.Intersections[5].Type);
            Assert.AreEqual(EdgeState.Filled, game.CurrentMesh.Edges[0].State);
            Assert.AreEqual(EdgeState.Excluded, game.CurrentMesh.Edges[1].State);

            // 3. Pearl Descriptor format
            string? desc = PuzzleFormatHelper.ExportPearlDescriptor(mesh);
            Assert.IsNotNull(desc);
            Assert.StartsWith("4x3:", desc);
            Assert.IsTrue(game.ImportPuzzle(desc, out string descFmt));
            Assert.AreEqual("Pearl Descriptor", descFmt);
            Assert.AreEqual(IntersType.Twist, game.CurrentMesh!.Intersections[0].Type);
            Assert.AreEqual(IntersType.Turn, game.CurrentMesh.Intersections[5].Type);

            // 4. Share Link
            string compact = PuzzleCodec.Encode(mesh, 3, 2);
            string fullUrl = $"https://www.themissingdocs.net/twistnturn/#p={compact}";
            Assert.IsTrue(game.ImportPuzzle(fullUrl, out string shareFmt));
            Assert.AreEqual("Share Link", shareFmt);
            Assert.AreEqual("3x2", game.SizeText);

            // 5. Rejects invalid text & rejects square formats for hexagonal grids
            Assert.IsFalse(game.ImportPuzzle("not a valid puzzle format!!!", out _));
            Assert.IsFalse(PuzzleFormatHelper.CanExportSquareFormats(new Mesh(5, 5, MeshType.Hexagonal2)));
            Assert.IsFalse(PuzzleFormatHelper.CanExportSquareFormats(new Mesh(3, 3, MeshType.Hexagonal4)));
        }

        [TestMethod]
        public void HandleEdgeGroupingClick_LinksSameAndOppositeParity_AndSupportsUndo()
        {
            var game = new GameService();
            game.LoadPuzzleText(PuzzleCodec.Encode(CreateSquare2x2CornerBlackPuzzle(), 2, 2));
            game.ShowEdgeColors = true;
            Assert.IsNotNull(game.CurrentMesh);

            // First click sets pending edge; clicking same edge cancels pending state
            Assert.IsTrue(game.HandleEdgeGroupingClick(0, sameLoop: true));
            Assert.AreEqual(0, game.PendingGroupEdgeIndex);
            Assert.IsTrue(game.HandleEdgeGroupingClick(0, sameLoop: true));
            Assert.IsNull(game.PendingGroupEdgeIndex);

            // Link edge 0 and edge 1 as same loop parity
            Assert.IsTrue(game.HandleEdgeGroupingClick(0, sameLoop: true));
            Assert.IsTrue(game.HandleEdgeGroupingClick(1, sameLoop: true));
            Assert.IsNull(game.PendingGroupEdgeIndex);
            Assert.AreNotEqual(0, game.CurrentMesh.Edges[0].Color);
            Assert.AreEqual(game.CurrentMesh.Edges[0].Color, game.CurrentMesh.Edges[1].Color);

            // Link edge 1 and edge 2 as opposite parity
            Assert.IsTrue(game.HandleEdgeGroupingClick(1, sameLoop: false));
            Assert.IsTrue(game.HandleEdgeGroupingClick(2, sameLoop: false));
            Assert.AreEqual(-game.CurrentMesh.Edges[1].Color, game.CurrentMesh.Edges[2].Color);

            // Contradictory link (0 and 2 are opposite, try forcing them same) should fail cleanly
            Assert.IsTrue(game.HandleEdgeGroupingClick(0, sameLoop: true));
            Assert.IsFalse(game.HandleEdgeGroupingClick(2, sameLoop: true));

            // Undo reverses the opposite-parity link cleanly
            game.Undo();
            Assert.AreEqual(0, game.CurrentMesh.Edges[2].Color);

            // Disabling ShowEdgeColors clears pending selection and round-trips in preferences
            Assert.IsTrue(game.HandleEdgeGroupingClick(0, sameLoop: true));
            Assert.AreEqual(0, game.PendingGroupEdgeIndex);
            var exported = game.ExportPreferences();
            Assert.IsTrue(exported.ShowEdgeColors);
            game.ShowEdgeColors = false;
            Assert.IsNull(game.PendingGroupEdgeIndex);
        }

        [TestMethod]
        public void AnalyzeCurrentPuzzle_ReportsUniqueMultipleAndUnsolvableCorrectly()
        {
            var game = new GameService();

            // 1. Unique solvable Masyu puzzle (4 Black corner pearls on 2x2)
            game.LoadPuzzleText(PuzzleCodec.Encode(CreateSquare2x2CornerBlackPuzzle(), 2, 2));
            var unique = game.AnalyzeCurrentPuzzle();
            Assert.IsNotNull(unique);
            Assert.AreEqual("Unique Solution", unique.SolveStatus);
            Assert.AreEqual(0, unique.BranchDepth);
            Assert.AreEqual(0, unique.WhitePearls);
            Assert.AreEqual(4, unique.BlackPearls);
            Assert.AreEqual(9, unique.TotalIntersections);
            Assert.StartsWith("Direct Deductions (S0)", unique.TechniqueTier);

            // 2. Empty grid -> Multiple Solutions
            var blank = new Mesh(2, 2, MeshType.Square);
            game.LoadPuzzleText(PuzzleCodec.Encode(blank, 2, 2));
            Assert.IsNull(game.CurrentAnalysis, "Loading new puzzle must clear prior analysis");
            var multiple = game.AnalyzeCurrentPuzzle();
            Assert.IsNotNull(multiple);
            Assert.AreEqual("Multiple Solutions", multiple.SolveStatus);

            // 3. Impossible Masyu clues -> No Solution (2x2 corner with White pearl cannot go straight through corner)
            var impossible = new Mesh(2, 2, MeshType.Square);
            impossible.Intersections[0].Type = IntersType.Twist; // Corner (0,0) has no straight edges!
            game.LoadPuzzleText(PuzzleCodec.Encode(impossible, 2, 2));
            var unsolvable = game.AnalyzeCurrentPuzzle();
            Assert.IsNotNull(unsolvable);
            Assert.AreEqual("No Solution", unsolvable.SolveStatus);
        }
    }
}
