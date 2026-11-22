using System.Text.Json;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using TwistNTurn;
using TwistNTurn.Web.Services;

namespace TwistNTurn.Tests
{
    [TestClass]
    public class Phase2SettingsTests
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
        public void PuzzleSetEdgeStateAction_SupportsAutoMoveAndEquality()
        {
            var mesh = CreateSquare2x2CornerBlackPuzzle();
            var a1 = new PuzzleSetEdgeStateAction(mesh, 0, EdgeState.Filled, disallowFalseMove: true, autoMove: 4);
            var a2 = new PuzzleSetEdgeStateAction(mesh, 0, EdgeState.Filled, disallowFalseMove: true, autoMove: 4);
            var a3 = new PuzzleSetEdgeStateAction(mesh, 0, EdgeState.Filled, disallowFalseMove: true, autoMove: 0);

            Assert.IsTrue(a1.Equals(a2));
            Assert.IsFalse(a1.Equals(a3));

            Assert.IsTrue(a1.Perform());
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);
            a1.Unperform();
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[0].State);
        }

        [TestMethod]
        public void TwistJsonContext_RoundTripsUserPreferencesAndSavedGameState()
        {
            var original = new UserPreferences
            {
                DisallowFalseMove = true,
                ShowCellColors = true,
                ShowEdgeColors = true,
                AutoMoveLevel = 3,
                AutoMoveConsiderMultipleLoops = false,
                AutoMoveUseColoring = true,
                AutoStart = true,
                WheelMode = MouseWheelMode.UndoRedo,
                ShowClock = false,
                GenerateLengthFraction = 0.65,
                GenerateBoringFraction = 0.05,
                SolverModeOverride = GeneratorSolverMode.Recursive
            };

            string json = JsonSerializer.Serialize(original, TwistJsonContext.Default.UserPreferences);
            var restored = JsonSerializer.Deserialize(json, TwistJsonContext.Default.UserPreferences);
            Assert.IsNotNull(restored);
            Assert.AreEqual(original.AutoMoveLevel, restored.AutoMoveLevel);
            Assert.AreEqual(original.AutoMoveUseColoring, restored.AutoMoveUseColoring);
            Assert.AreEqual(original.WheelMode, restored.WheelMode);
            Assert.AreEqual(original.ShowEdgeColors, restored.ShowEdgeColors);
            Assert.AreEqual(original.SolverModeOverride, restored.SolverModeOverride);

            using var game = new GameService();
            game.ApplyPreferences(restored);
            Assert.IsTrue(game.AutoStart);
            Assert.IsFalse(game.ShowClock);

            game.ResetPreferencesToDefaults();
            Assert.AreEqual(0, game.AutoMoveLevel);
            Assert.IsTrue(game.ShowClock);
            Assert.IsFalse(game.ShowEdgeColors);
        }

        [TestMethod]
        public void AutoStart_AndGeneratorOverrides_ApplyCleanly()
        {
            using var game = new GameService();
            game.AutoStart = true;
            game.LoadPuzzleText(PuzzleCodec.Encode(CreateSquare2x2CornerBlackPuzzle(), 2, 2));

            // Reset triggers ApplyAutoStartIfEnabled
            game.ResetCurrentPuzzle();
            int nonEmpty = 0;
            foreach (var e in game.CurrentMesh!.Edges)
            {
                if (e.State != EdgeState.Empty) nonEmpty++;
            }
            Assert.IsGreaterThan(0, nonEmpty, "Auto-start should fill forced corner pearl edges.");

            game.GenerateLengthFraction = 0.7;
            game.GenerateBoringFraction = 0.08;
            game.SolverModeOverride = GeneratorSolverMode.Iterative;
            var target = new Mesh(4, 4, MeshType.Square);
            game.ApplyGeneratorOverrides(target);
            Assert.AreEqual(0.7, target.GenerateLengthFraction, 1e-6);
            Assert.AreEqual(0.08, target.GenerateBoringFraction, 1e-6);
            Assert.AreEqual(SolverMethod.Iterative, target.SolverMethod);
        }
    }
}
