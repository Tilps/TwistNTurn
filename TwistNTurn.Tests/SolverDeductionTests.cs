using System;
using System.Threading;
using System.Threading.Tasks;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using TwistNTurn;

namespace TwistNTurn.Tests
{
    [TestClass]
    public class SolverDeductionTests
    {
        [TestMethod]
        public async Task Generator_ProducesSolvablePuzzle_OnSmallGrid()
        {
            var mesh = new Mesh(3, 3, MeshType.Square);
            using var cts = new CancellationTokenSource(TimeSpan.FromSeconds(30));

            await mesh.GenerateAsync(progress: null, cts.Token);

            // Verify some pearls were placed
            int clueCount = 0;
            for (int i = 0; i < mesh.Intersections.Count; i++)
            {
                if (mesh.Intersections[i].Type != IntersType.Unknown)
                    clueCount++;
            }
            Assert.IsTrue(clueCount > 0, "Generator should produce at least one clue pearl");

            // Verify the generated puzzle solves to Solved
            mesh.Clear(); // reset edges to empty
            var solveState = mesh.TrySolve();
            Assert.AreEqual(SolveState.Solved, solveState, "Generated puzzle must be uniquely solvable");
        }

        [TestMethod]
        public void Solver_ContradictionPuzzle_ReturnsNoSolutions()
        {
            var mesh = new Mesh(2, 2, MeshType.Square);
            // Place conflicting clues: on a 2x2 grid, all intersections turn on 90 degrees if they touch filled edges.
            // A Turn (black pearl) demands both adjacent edges go straight. On a 2x2 grid, no straight path of length 2 exists.
            // Placing Turn pearls on corners will make it impossible to solve.
            mesh.Intersections[0].Type = IntersType.Turn;
            mesh.Intersections[1].Type = IntersType.Turn;

            var solveState = mesh.TrySolve();
            Assert.AreEqual(SolveState.NoSolutions, solveState);
        }
    }
}

