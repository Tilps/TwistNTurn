using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using TwistNTurn;

namespace TwistNTurn.Tests
{
    [TestClass]
    public class ActionHistoryTests
    {
        [TestMethod]
        public void PuzzleEdgeAction_NormalCycle_TransitionsAndUndosCorrectly()
        {
            var mesh = new Mesh(3, 3, MeshType.Square);
            var tree = new UndoTree();

            Assert.AreEqual(EdgeState.Empty, mesh.Edges[0].State);

            // Empty -> Filled
            var a1 = new PuzzleEdgeAction(mesh, 0, isAlternativeCycle: false);
            Assert.IsTrue(tree.Do(a1));
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);

            // Filled -> Excluded
            var a2 = new PuzzleEdgeAction(mesh, 0, isAlternativeCycle: false);
            Assert.IsTrue(tree.Do(a2));
            Assert.AreEqual(EdgeState.Excluded, mesh.Edges[0].State);

            // Excluded -> Empty
            var a3 = new PuzzleEdgeAction(mesh, 0, isAlternativeCycle: false);
            Assert.IsTrue(tree.Do(a3));
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[0].State);

            // Undo all 3 steps
            Assert.IsTrue(tree.Undo());
            Assert.AreEqual(EdgeState.Excluded, mesh.Edges[0].State);

            Assert.IsTrue(tree.Undo());
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);

            Assert.IsTrue(tree.Undo());
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[0].State);

            // Redo all 3 steps
            Assert.IsTrue(tree.Redo());
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);

            Assert.IsTrue(tree.Redo());
            Assert.AreEqual(EdgeState.Excluded, mesh.Edges[0].State);

            Assert.IsTrue(tree.Redo());
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[0].State);
        }

        [TestMethod]
        public void PuzzleEdgeAction_AlternativeCycle_TransitionsCorrectly()
        {
            var mesh = new Mesh(3, 3, MeshType.Square);
            var tree = new UndoTree();

            // Alternative cycle (right-click / alt mode): Empty -> Excluded -> Filled -> Empty
            var a1 = new PuzzleEdgeAction(mesh, 0, isAlternativeCycle: true);
            Assert.IsTrue(tree.Do(a1));
            Assert.AreEqual(EdgeState.Excluded, mesh.Edges[0].State);

            var a2 = new PuzzleEdgeAction(mesh, 0, isAlternativeCycle: true);
            Assert.IsTrue(tree.Do(a2));
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);

            var a3 = new PuzzleEdgeAction(mesh, 0, isAlternativeCycle: true);
            Assert.IsTrue(tree.Do(a3));
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[0].State);

            // Undo step
            Assert.IsTrue(tree.Undo());
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);
        }

        [TestMethod]
        public void PuzzleSetEdgeStateAction_SetsDirectStateAndUndos()
        {
            var mesh = new Mesh(3, 3, MeshType.Square);
            var tree = new UndoTree();

            var action = new PuzzleSetEdgeStateAction(mesh, 0, EdgeState.Filled);
            Assert.IsTrue(tree.Do(action));
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);

            // Performing identical state should be a no-op
            var noop = new PuzzleSetEdgeStateAction(mesh, 0, EdgeState.Filled);
            Assert.IsFalse(tree.Do(noop));

            // Undo
            Assert.IsTrue(tree.Undo());
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[0].State);
        }

        [TestMethod]
        public void UndoTree_Checkpoints_FixAndRevertWorkCorrectly()
        {
            var mesh = new Mesh(3, 3, MeshType.Square);
            var tree = new UndoTree();

            // Perform two moves
            tree.Do(new PuzzleEdgeAction(mesh, 0, false)); // edge 0 filled
            tree.Do(new PuzzleEdgeAction(mesh, 1, false)); // edge 1 filled

            Assert.IsFalse(tree.IsMarked);
            Assert.IsFalse(tree.CanRevert);

            // Fix checkpoint
            tree.Mark();
            Assert.IsTrue(tree.IsMarked);
            Assert.IsFalse(tree.CanRevert); // no moves made since mark

            // Make two more moves
            tree.Do(new PuzzleEdgeAction(mesh, 2, false)); // edge 2 filled
            tree.Do(new PuzzleEdgeAction(mesh, 3, true));  // edge 3 excluded

            Assert.IsTrue(tree.CanRevert);
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[2].State);
            Assert.AreEqual(EdgeState.Excluded, mesh.Edges[3].State);

            // Revert to mark
            tree.RevertToMark();
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[2].State);
            Assert.AreEqual(EdgeState.Empty, mesh.Edges[3].State);
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[0].State);
            Assert.AreEqual(EdgeState.Filled, mesh.Edges[1].State);

            // Clear mark (Unfix)
            tree.ClearMark();
            Assert.IsFalse(tree.IsMarked);
            Assert.IsFalse(tree.CanRevert);
        }

        [TestMethod]
        public void PuzzleCellColorAction_CyclesColorsAndUndos()
        {
            var mesh = new Mesh(3, 3, MeshType.Square);
            var tree = new UndoTree();

            Assert.AreEqual(0, mesh.Cells[0].Color);

            // Set cell color from 0 to 1
            tree.Do(new PuzzleCellColorAction(mesh, 0, reverse: false));
            Assert.AreEqual(1, mesh.Cells[0].Color);

            // Set cell color from 1 to -1
            tree.Do(new PuzzleCellColorAction(mesh, 0, reverse: false));
            Assert.AreEqual(-1, mesh.Cells[0].Color);

            // Undo to 1
            Assert.IsTrue(tree.Undo());
            Assert.AreEqual(1, mesh.Cells[0].Color);

            // Undo to 0
            Assert.IsTrue(tree.Undo());
            Assert.AreEqual(0, mesh.Cells[0].Color);
        }
    }
}
