using Microsoft.VisualStudio.TestTools.UnitTesting;
using TwistNTurn;
using TwistNTurn.Web.Services;

namespace TwistNTurn.Tests
{
    [TestClass]
    public class Phase3VisualSolverTests
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
        public void StartVisualSolve_RecordsFrames_Navigates_AndAppliesSolution()
        {
            var game = new GameService();
            game.LoadPuzzleText(PuzzleCodec.Encode(CreateSquare2x2CornerBlackPuzzle(), 2, 2));

            Assert.IsTrue(game.StartVisualSolve());
            Assert.IsTrue(game.IsVisualSolving);
            Assert.IsGreaterThanOrEqualTo(2, game.VisualSolveTotalSteps);
            Assert.AreEqual(0, game.VisualSolveStepIndex);

            // Step navigation & clamping
            Assert.IsFalse(game.StepVisualSolveBackward());
            Assert.IsTrue(game.StepVisualSolveForward());
            Assert.AreEqual(1, game.VisualSolveStepIndex);

            game.SetVisualSolveStep(9999);
            Assert.AreEqual(game.VisualSolveTotalSteps - 1, game.VisualSolveStepIndex);
            Assert.AreEqual("Solution complete!", game.CurrentVisualSolveStep!.Status);
            Assert.IsFalse(game.StepVisualSolveForward());

            // Applying solution solves puzzle and exits visual solve mode
            Assert.IsTrue(game.ApplyVisualSolveSolution());
            Assert.IsFalse(game.IsVisualSolving);
            Assert.IsTrue(game.IsSolved);
            Assert.IsTrue(game.WasAutoSolved);

            // Undo reverses the applied solution cleanly
            game.Undo();
            Assert.IsFalse(game.IsSolved);
        }

        [TestMethod]
        public void DismissVisualSolve_LeavesOriginalBoardUntouched()
        {
            var game = new GameService();
            game.LoadPuzzleText(PuzzleCodec.Encode(CreateSquare2x2CornerBlackPuzzle(), 2, 2));

            Assert.IsTrue(game.StartVisualSolve());
            game.SetVisualSolveStep(game.VisualSolveTotalSteps - 1);
            game.DismissVisualSolve();

            Assert.IsFalse(game.IsVisualSolving);
            Assert.IsFalse(game.IsSolved);
            Assert.AreEqual(EdgeState.Empty, game.CurrentMesh!.Edges[0].State);
        }
    }
}
