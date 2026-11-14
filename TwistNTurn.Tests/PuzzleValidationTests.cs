using System;
using System.Threading.Tasks;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using TwistNTurn;

namespace TwistNTurn.Tests
{
    [TestClass]
    public class PuzzleValidationTests
    {
        [TestMethod]
        public void CheckIsSolved_EmptyBoard_ReturnsFalse()
        {
            var mesh = new Mesh(4, 4, MeshType.Square);
            Assert.IsFalse(PuzzleHelper.CheckIsSolved(mesh));
        }

        [TestMethod]
        public void CheckIsSolved_Simple2x2Loop_ValidatesCorrectly()
        {
            // Create a 2x2 square mesh
            var mesh = new Mesh(2, 2, MeshType.Square);

            // A 2x2 mesh of squares has (2+1)*(2+1) = 9 intersections, or 4 cells.
            // Let's form a single outer loop around cell 0 (a single unit square)
            // Let's find the 4 edges of cell 0
            var cell0 = mesh.Cells[0];
            Assert.AreEqual(4, cell0.Edges.Count);

            foreach (int edgeIndex in cell0.Edges)
            {
                mesh.Edges[edgeIndex].State = EdgeState.Filled;
            }

            // Exclude other edges
            for (int i = 0; i < mesh.Edges.Count; i++)
            {
                if (mesh.Edges[i].State == EdgeState.Empty)
                    mesh.Edges[i].State = EdgeState.Excluded;
            }

            // Without clues, any closed cycle of degree 2 is solved!
            Assert.IsTrue(PuzzleHelper.CheckIsSolved(mesh));
        }

        [TestMethod]
        public void DetermineType_TurnPearl_RulesValidation()
        {
            // Build a 3x3 square mesh (4x4 intersections)
            var mesh = new Mesh(3, 3, MeshType.Square);

            // In a loop of a single unit square:
            // All 4 vertices turn 90 degrees.
            // For a Turn (black pearl), both adjacent vertices MUST go straight.
            // Because the corners of a 1x1 square turn immediately on both sides,
            // placing a Turn (black) on a 1x1 square is INVALID.
            var cell0 = mesh.Cells[0];
            foreach (int edgeIndex in cell0.Edges)
            {
                mesh.Edges[edgeIndex].State = EdgeState.Filled;
            }

            // Cell 0 corners turn on both adjacent vertices
            int cornerInters = cell0.Intersections[0];
            mesh.Intersections[cornerInters].Type = IntersType.Turn;
            mesh.RecalculateCounts();

            // Turn rule is NOT satisfied because adjacent vertices turn instead of going straight
            Assert.IsFalse(mesh.TypeSatisfied(cornerInters));
            Assert.IsFalse(PuzzleHelper.CheckIsSolved(mesh));
        }

        [TestMethod]
        public void DetermineType_TwistPearl_RulesValidation()
        {
            // Build a 3x3 square mesh (4x4 intersections)
            var mesh = new Mesh(3, 3, MeshType.Square);

            // Let's create a 2x1 rectangle loop:
            // Intersections:
            // (0,0) -> (1,0) -> (2,0) -> (2,1) -> (1,1) -> (0,1) -> (0,0)
            // Vertices (1,0) and (1,1) are straight, with adjacent vertices (0,0)/(2,0) and (0,1)/(2,1) turning!
            // Therefore, (1,0) is an ideal Twist (white pearl) candidate!
            // Let's fill the edges of cell 0 and cell 1 (adjacent horizontally)
            var cell0 = mesh.Cells[0];
            var cell1 = mesh.Cells[1];

            // Fill all outer edges of the 2x1 rectangle
            foreach (int e in cell0.Edges)
                mesh.Edges[e].State = EdgeState.Filled;
            foreach (int e in cell1.Edges)
                mesh.Edges[e].State = EdgeState.Filled;

            // Remove the internal shared dividing edge
            int sharedEdge = -1;
            foreach (int e in cell0.Edges)
            {
                if (cell1.Edges.Contains(e))
                {
                    sharedEdge = e;
                    break;
                }
            }
            Assert.AreNotEqual(-1, sharedEdge);
            mesh.Edges[sharedEdge].State = EdgeState.Excluded;

            // Find an intersection along the straight edge
            int straightInters = -1;
            foreach (int intersIndex in mesh.Edges[sharedEdge].Intersections)
            {
                // Count filled edges at this intersection
                int filled = 0;
                foreach (int e in mesh.Intersections[intersIndex].Edges)
                {
                    if (mesh.Edges[e].State == EdgeState.Filled) filled++;
                }
                if (filled == 2)
                {
                    straightInters = intersIndex;
                    break;
                }
            }
            Assert.AreNotEqual(-1, straightInters);

            // Make straightInters a Twist (white pearl)
            mesh.Intersections[straightInters].Type = IntersType.Twist;
            mesh.RecalculateCounts();

            // Verify Twist is satisfied because it goes straight and adjacent corners turn!
            Assert.IsTrue(mesh.TypeSatisfied(straightInters), "Twist pearl should be satisfied");
            Assert.AreEqual(IntersType.Twist, mesh.DetermineType(straightInters));
            Assert.IsTrue(PuzzleHelper.CheckIsSolved(mesh));
        }

        [TestMethod]
        public void CheckIsSolved_DisconnectedTwoLoops_ReturnsFalse()
        {
            // Build a 4x4 square mesh (enough room for two separate 1x1 loops)
            var mesh = new Mesh(4, 4, MeshType.Square);

            // Cell 0 (top-left)
            foreach (int e in mesh.Cells[0].Edges)
                mesh.Edges[e].State = EdgeState.Filled;

            // Cell at bottom-right (no shared edges or intersections with cell 0)
            int lastCell = mesh.Cells.Count - 1;
            foreach (int e in mesh.Cells[lastCell].Edges)
                mesh.Edges[e].State = EdgeState.Filled;

            // Even though all filled intersections have degree 2, there are 2 disconnected cycles
            Assert.IsFalse(PuzzleHelper.CheckIsSolved(mesh),
                "Two disconnected loops should NOT be considered solved");
        }

        [TestMethod]
        public async Task GenerateAsync_CancelledImmediately_AbortsCleanly()
        {
            var mesh = PuzzleHelper.MakeMesh(3, 3, MeshType.Square, 0);
            using var cts = new System.Threading.CancellationTokenSource();
            cts.Cancel();
            await mesh.GenerateAsync(cancellationToken: cts.Token);
            Assert.IsTrue(cts.IsCancellationRequested);
        }

        [TestMethod]
        public async Task GenerateAsync_AbortedViaAbortPrune_AbortsCleanly()
        {
            var mesh = PuzzleHelper.MakeMesh(6, 6, MeshType.Square, 0);
            var task = mesh.GenerateAsync();
            mesh.AbortPrune = true;
            await task;
            Assert.IsTrue(mesh.AbortPrune);
        }

        [TestMethod]
        [DataRow(MeshType.Square, 3, 3)]
        [DataRow(MeshType.Hexagonal2, 5, 5)]
        [DataRow(MeshType.Hexagonal4, 3, 3)]
        public async Task GenerateAsync_ProducesValidPuzzle(MeshType type, int width, int height)
        {
            var mesh = PuzzleHelper.MakeMesh(width, height, type, 0);
            await mesh.GenerateAsync();
            Assert.AreEqual(SolveState.Solved, mesh.TrySolve());
        }
    }
}

