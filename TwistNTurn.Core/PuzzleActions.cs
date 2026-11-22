using System;
using System.Collections.Generic;

namespace TwistNTurn
{
    public class PuzzleEdgeAction : IAction
    {
        public PuzzleEdgeAction(Mesh mesh, int edgeIndex, bool isAlternativeCycle = false, int autoMove = 0, bool disallowFalseMove = false)
        {
            this.mesh = mesh;
            this.edgeIndex = edgeIndex;
            this.isAlternativeCycle = isAlternativeCycle;
            this.autoMove = autoMove;
            this.disallowFalseMove = disallowFalseMove;
        }

        private readonly Mesh mesh;
        private readonly int edgeIndex;
        private readonly bool isAlternativeCycle;
        private readonly int autoMove;
        private readonly bool disallowFalseMove;
        private List<IAction>? actionsPerformed;
        private bool successful;

        public bool Successful => successful;
        public string Name => $"Edge {edgeIndex} Toggle";

        internal static EdgeState NextState(EdgeState current, bool reverse) => (current, reverse) switch
        {
            (EdgeState.Empty, false) => EdgeState.Filled,
            (EdgeState.Filled, false) => EdgeState.Excluded,
            (EdgeState.Empty, true) => EdgeState.Excluded,
            (EdgeState.Excluded, true) => EdgeState.Filled,
            _ => EdgeState.Empty
        };

        internal static bool TransitionEdge(
            Mesh mesh,
            int edgeIndex,
            EdgeState targetState,
            int autoMove,
            bool disallowFalseMove,
            out List<IAction> actionsPerformed,
            out bool successful)
        {
            successful = true;
            actionsPerformed = [];
            Edge closest = mesh.Edges[edgeIndex];

            if (closest.State != EdgeState.Empty)
            {
                var unsetAction = new UnsetAction(mesh, edgeIndex);
                bool res = unsetAction.Perform();
                if (!res || !unsetAction.Successful)
                {
                    successful = false;
                    if (res)
                        actionsPerformed.Add(unsetAction);
                    if (disallowFalseMove)
                    {
                        if (actionsPerformed.Count > 0)
                        {
                            mesh.Unperform(actionsPerformed);
                            actionsPerformed.Clear();
                        }
                        return false;
                    }
                }
                else
                {
                    actionsPerformed.Add(unsetAction);
                }
            }

            if (targetState != EdgeState.Empty)
            {
                bool res = mesh.Perform(edgeIndex, targetState, actionsPerformed, autoMove);
                if (!res)
                {
                    successful = false;
                    if (disallowFalseMove)
                    {
                        if (actionsPerformed.Count > 0)
                        {
                            mesh.Unperform(actionsPerformed);
                            actionsPerformed.Clear();
                        }
                        return false;
                    }
                }
            }

            return true;
        }

        public bool Perform()
        {
            EdgeState toggled = NextState(mesh.Edges[edgeIndex].State, isAlternativeCycle);
            return TransitionEdge(mesh, edgeIndex, toggled, autoMove, disallowFalseMove, out actionsPerformed, out successful);
        }

        public void Unperform()
        {
            if (actionsPerformed is { Count: > 0 })
            {
                mesh.Unperform(actionsPerformed);
                actionsPerformed.Clear();
            }
        }

        public bool Equals(IAction? other)
        {
            return other is PuzzleEdgeAction o
                && o.mesh == mesh
                && o.edgeIndex == edgeIndex
                && o.isAlternativeCycle == isAlternativeCycle
                && o.autoMove == autoMove
                && o.disallowFalseMove == disallowFalseMove;
        }

        public override bool Equals(object? obj) => obj is IAction action && Equals(action);

        public override int GetHashCode() => System.HashCode.Combine(mesh, edgeIndex, isAlternativeCycle, autoMove, disallowFalseMove);
    }

    public class PuzzleCellColorAction : IAction
    {
        public PuzzleCellColorAction(Mesh mesh, int cellIndex, bool reverse = false)
        {
            this.mesh = mesh;
            this.cellIndex = cellIndex;
            this.reverse = reverse;
        }

        private readonly Mesh mesh;
        private readonly int cellIndex;
        private readonly bool reverse;
        private List<IAction>? actionsPerformed;
        private bool successful;

        public bool Successful => successful;
        public string Name => $"Cell {cellIndex} Color";

        private static int Toggle(int color, bool reverse)
        {
            if (!reverse)
            {
                if (color == 0) return 1;
                if (color == 1) return -1;
                return 0;
            }
            else
            {
                if (color == 0) return -1;
                if (color == -1) return 1;
                return 0;
            }
        }

        public bool Perform()
        {
            successful = false;
            Cell closest = mesh.Cells[cellIndex];
            int newColor = Toggle(closest.Color, reverse);
            actionsPerformed = [];

            if (closest.Color != 0)
            {
                var unsetAction = new CellColorClearAction(mesh, cellIndex);
                bool res = unsetAction.Perform();
                if (!res || !unsetAction.Successful)
                {
                    if (res)
                        actionsPerformed.Add(unsetAction);
                    Unperform();
                    return false;
                }
                actionsPerformed.Add(unsetAction);
            }

            if (newColor != 0)
            {
                var setAction = new CellColorJoinAction(mesh, cellIndex, -1, newColor == 1);
                bool res = setAction.Perform();
                if (!res || !setAction.Successful)
                {
                    if (res)
                        actionsPerformed.Add(setAction);
                    Unperform();
                    return false;
                }
                actionsPerformed.Add(setAction);
            }

            successful = true;
            return true;
        }

        public void Unperform()
        {
            if (actionsPerformed is { Count: > 0 })
            {
                mesh.Unperform(actionsPerformed);
                actionsPerformed.Clear();
            }
        }

        public bool Equals(IAction? other)
        {
            return other is PuzzleCellColorAction o && o.mesh == mesh && o.cellIndex == cellIndex && o.reverse == reverse;
        }

        public override bool Equals(object? obj) => obj is IAction action && Equals(action);

        public override int GetHashCode() => System.HashCode.Combine(mesh, cellIndex, reverse);
    }

    public class PuzzleSetEdgeStateAction : IAction
    {
        public PuzzleSetEdgeStateAction(Mesh mesh, int edgeIndex, EdgeState targetState, bool disallowFalseMove = false, int autoMove = 0)
        {
            this.mesh = mesh;
            this.edgeIndex = edgeIndex;
            this.targetState = targetState;
            this.disallowFalseMove = disallowFalseMove;
            this.autoMove = autoMove;
        }

        private readonly Mesh mesh;
        private readonly int edgeIndex;
        private readonly EdgeState targetState;
        private readonly bool disallowFalseMove;
        private readonly int autoMove;
        private List<IAction>? actionsPerformed;
        private bool successful;

        public bool Successful => successful;
        public string Name => $"Edge {edgeIndex} -> {targetState}";

        public bool Perform()
        {
            if (mesh.Edges[edgeIndex].State == targetState)
            {
                successful = true;
                return false;
            }
            return PuzzleEdgeAction.TransitionEdge(mesh, edgeIndex, targetState, autoMove, disallowFalseMove, out actionsPerformed, out successful);
        }

        public void Unperform()
        {
            if (actionsPerformed is { Count: > 0 })
            {
                mesh.Unperform(actionsPerformed);
                actionsPerformed.Clear();
            }
        }

        public bool Equals(IAction? other)
        {
            return other is PuzzleSetEdgeStateAction o
                && o.mesh == mesh
                && o.edgeIndex == edgeIndex
                && o.targetState == targetState
                && o.disallowFalseMove == disallowFalseMove
                && o.autoMove == autoMove;
        }

        public override bool Equals(object? obj) => obj is IAction action && Equals(action);

        public override int GetHashCode() => System.HashCode.Combine(mesh, edgeIndex, targetState, disallowFalseMove, autoMove);
    }

    public class PuzzleBatchSolveAction : IAction
    {
        private readonly Mesh mesh;
        private readonly EdgeState[] newStates;
        private readonly List<IAction> performed = [];

        public PuzzleBatchSolveAction(Mesh mesh, Mesh solvedMesh)
        {
            if (solvedMesh.Edges.Count != mesh.Edges.Count)
                throw new ArgumentException("Solved mesh edge count does not match target mesh.", nameof(solvedMesh));

            this.mesh = mesh;
            newStates = new EdgeState[mesh.Edges.Count];
            for (int i = 0; i < mesh.Edges.Count; i++)
            {
                newStates[i] = solvedMesh.Edges[i].State;
            }
        }

        public bool Successful => true;
        public string Name => "Solve Puzzle";

        public bool Perform()
        {
            performed.Clear();
            for (int i = 0; i < mesh.Edges.Count; i++)
            {
                if (mesh.Edges[i].State != EdgeState.Empty)
                {
                    var unset = new UnsetAction(mesh, i);
                    unset.Perform();
                    performed.Add(unset);
                }
            }
            for (int i = 0; i < mesh.Edges.Count; i++)
            {
                if (newStates[i] != EdgeState.Empty)
                {
                    var set = new SetAction(mesh, i, newStates[i]);
                    set.Perform();
                    performed.Add(set);
                }
            }
            return true;
        }

        public void Unperform()
        {
            for (int i = performed.Count - 1; i >= 0; i--)
            {
                performed[i].Unperform();
            }
            performed.Clear();
        }

        public bool Equals(IAction? other)
        {
            return other is PuzzleBatchSolveAction o
                && o.mesh == mesh
                && o.newStates.AsSpan().SequenceEqual(newStates);
        }

        public override bool Equals(object? obj) => obj is IAction action && Equals(action);

        public override int GetHashCode()
        {
            var hc = new System.HashCode();
            hc.Add(mesh);
            for (int i = 0; i < newStates.Length; i++)
                hc.Add(newStates[i]);
            return hc.ToHashCode();
        }
    }
}
