using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using TwistNTurn;

namespace TwistNTurn.Tests
{
    [TestClass]
    public class TopologyTests
    {
        [TestMethod]
        [DataRow(MeshType.Square, 4, 4)]
        [DataRow(MeshType.Square, 6, 6)]
        [DataRow(MeshType.SquareSymmetrical, 4, 4)]
        [DataRow(MeshType.SquareSymmetrical, 6, 6)]
        [DataRow(MeshType.Hexagonal3, 3, 3)]
        public void Topologies_ConstructGrid_InitializesValidGraph(MeshType type, int width, int height)
        {
            var mesh = new Mesh(width, height, type);

            Assert.IsNotNull(mesh.Edges);
            Assert.IsNotNull(mesh.Intersections);
            Assert.IsNotNull(mesh.Cells);

            Assert.IsTrue(mesh.Edges.Count > 0, $"Edges should exist for {type}");
            Assert.IsTrue(mesh.Intersections.Count > 0, $"Intersections should exist for {type}");
            Assert.IsTrue(mesh.Cells.Count > 0, $"Cells should exist for {type}");

            // Verify edges connect two valid intersections
            for (int i = 0; i < mesh.Edges.Count; i++)
            {
                var edge = mesh.Edges[i];
                Assert.AreEqual(2, edge.Intersections.Length, $"Edge {i} must have 2 intersections");
                Assert.IsTrue(edge.Intersections[0] >= 0 && edge.Intersections[0] < mesh.Intersections.Count,
                    $"Edge {i} intersection 0 is out of bounds");
                Assert.IsTrue(edge.Intersections[1] >= 0 && edge.Intersections[1] < mesh.Intersections.Count,
                    $"Edge {i} intersection 1 is out of bounds");
                Assert.AreNotEqual(edge.Intersections[0], edge.Intersections[1],
                    $"Edge {i} must connect distinct intersections");
            }

            // Verify opposite edges are populated for all intersections
            for (int i = 0; i < mesh.Intersections.Count; i++)
            {
                var inters = mesh.Intersections[i];
                Assert.AreEqual(inters.Edges.Count, inters.Opposite.Count,
                    $"Intersection {i} Opposites list count must match Edges count");

                for (int j = 0; j < inters.Edges.Count; j++)
                {
                    int edgeIndex = inters.Edges[j];
                    Assert.IsTrue(edgeIndex >= 0 && edgeIndex < mesh.Edges.Count,
                        $"Intersection {i} has invalid edge index {edgeIndex}");

                    int oppIndex = inters.Opposite[j];
                    if (oppIndex != -1)
                    {
                        Assert.IsTrue(oppIndex >= 0 && oppIndex < mesh.Edges.Count,
                            $"Intersection {i} has invalid opposite edge {oppIndex}");
                        Assert.AreNotEqual(edgeIndex, oppIndex,
                            $"Edge cannot be its own opposite at intersection {i}");
                    }
                }
            }

            // Verify cell definitions
            for (int i = 0; i < mesh.Cells.Count; i++)
            {
                var cell = mesh.Cells[i];
                Assert.IsTrue(cell.Edges.Count >= 3, $"Cell {i} must have at least 3 edges");
                Assert.AreEqual(cell.Edges.Count, cell.Intersections.Count,
                    $"Cell {i} edges and intersections count must match");

                for (int j = 0; j < cell.Edges.Count; j++)
                {
                    int edgeIndex = cell.Edges[j];
                    var edge = mesh.Edges[edgeIndex];
                    Assert.IsTrue(edge.Cells.Contains(i),
                        $"Edge {edgeIndex} must reference Cell {i}");
                }
            }
        }

        [TestMethod]
        public void SquareGrid_Opposites_FormStraightLines()
        {
            // In a 4x4 square grid, internal intersections have degree 4 and 2 pairs of opposites
            var mesh = new Mesh(4, 4, MeshType.Square);

            int internalIntersWithFourEdges = 0;
            for (int i = 0; i < mesh.Intersections.Count; i++)
            {
                var inters = mesh.Intersections[i];
                if (inters.Edges.Count == 4)
                {
                    internalIntersWithFourEdges++;
                    // Each of the 4 edges should have an opposite
                    for (int j = 0; j < inters.Edges.Count; j++)
                    {
                        Assert.AreNotEqual(-1, inters.Opposite[j],
                            $"Internal intersection {i} edge {inters.Edges[j]} should have an opposite");
                    }
                }
            }

            Assert.IsTrue(internalIntersWithFourEdges > 0, "Must have internal intersections in 4x4 square mesh");
        }
    }
}

