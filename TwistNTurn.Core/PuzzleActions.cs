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

        private static EdgeState NextState(EdgeState current, bool reverse)
        {
            if (!reverse)
            {
                return current switch
                {
                    EdgeState.Empty => EdgeState.Filled,
                    EdgeState.Filled => EdgeState.Excluded,
                    EdgeState.Excluded => EdgeState.Empty,
                    _ => EdgeState.Empty
                };
            }
            else
            {
                return current switch
                {
                    EdgeState.Empty => EdgeState.Excluded,
                    EdgeState.Excluded => EdgeState.Filled,
                    EdgeState.Filled => EdgeState.Empty,
                    _ => EdgeState.Empty
                };
            }
        }

        public bool Perform()
        {
            successful = true;
            Edge closest = mesh.Edges[edgeIndex];
            EdgeState toggled = NextState(closest.State, isAlternativeCycle);
            actionsPerformed = new List<IAction>();

            if (closest.State != EdgeState.Empty)
            {
                IAction unsetAction = new UnsetAction(mesh, edgeIndex);
                bool res = unsetAction.Perform();
                if ((!res || !unsetAction.Successful) && disallowFalseMove)
                {
                    if (res && !unsetAction.Successful)
                        actionsPerformed.Add(unsetAction);
                    Unperform();
                    return false;
                }
                else if (!res || !unsetAction.Successful)
                {
                    successful = false;
                }
                actionsPerformed.Add(unsetAction);
            }

            if (toggled != EdgeState.Empty)
            {
                bool res = mesh.Perform(edgeIndex, toggled, actionsPerformed, autoMove);
                if (!res && disallowFalseMove)
                {
                    Unperform();
                    return false;
                }
                else if (!res)
                {
                    successful = false;
                }
            }

            return true;
        }

        public void Unperform()
        {
            if (actionsPerformed != null && actionsPerformed.Count > 0)
            {
                mesh.Unperform(actionsPerformed);
            }
        }

        public bool Equals(IAction? other)
        {
            if (other is not PuzzleEdgeAction o) return false;
            return o.mesh == mesh && o.edgeIndex == edgeIndex && o.isAlternativeCycle == isAlternativeCycle;
        }
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
                if (color == 1) return -1;
                if (color == -1) return 0;
                return 1;
            }
            else
            {
                if (color == 1) return 0;
                if (color == -1) return 1;
                return -1;
            }
        }

        public bool Perform()
        {
            successful = true;
            Cell closest = mesh.Cells[cellIndex];
            int newColor = Toggle(closest.Color, reverse);
            actionsPerformed = new List<IAction>();

            if (closest.Color != 0)
            {
                IAction unsetAction = new CellColorClearAction(mesh, cellIndex);
                bool res = unsetAction.Perform();
                if (!res || !unsetAction.Successful)
                {
                    if (res && !unsetAction.Successful)
                        actionsPerformed.Add(unsetAction);
                    Unperform();
                    return false;
                }
                actionsPerformed.Add(unsetAction);
            }

            if (newColor != 0)
            {
                IAction setAction = new CellColorJoinAction(mesh, cellIndex, -1, newColor == 1);
                bool res = setAction.Perform();
                if (!res || !setAction.Successful)
                {
                    if (res && !setAction.Successful)
                        actionsPerformed.Add(setAction);
                    Unperform();
                    return false;
                }
                actionsPerformed.Add(setAction);
            }

            return true;
        }

        public void Unperform()
        {
            if (actionsPerformed != null && actionsPerformed.Count > 0)
            {
                mesh.Unperform(actionsPerformed);
            }
        }

        public bool Equals(IAction? other)
        {
            if (other is not PuzzleCellColorAction o) return false;
            return o.mesh == mesh && o.cellIndex == cellIndex && o.reverse == reverse;
        }
    }

    public class PuzzleSetEdgeStateAction : IAction
    {
        public PuzzleSetEdgeStateAction(Mesh mesh, int edgeIndex, EdgeState targetState, bool disallowFalseMove = false)
        {
            this.mesh = mesh;
            this.edgeIndex = edgeIndex;
            this.targetState = targetState;
            this.disallowFalseMove = disallowFalseMove;
        }

        private readonly Mesh mesh;
        private readonly int edgeIndex;
        private readonly EdgeState targetState;
        private readonly bool disallowFalseMove;
        private List<IAction>? actionsPerformed;
        private bool successful;

        public bool Successful => successful;
        public string Name => $"Edge {edgeIndex} -> {targetState}";

        public bool Perform()
        {
            successful = true;
            Edge closest = mesh.Edges[edgeIndex];
            if (closest.State == targetState) return false;

            actionsPerformed = new List<IAction>();

            if (closest.State != EdgeState.Empty)
            {
                IAction unsetAction = new UnsetAction(mesh, edgeIndex);
                bool res = unsetAction.Perform();
                if ((!res || !unsetAction.Successful) && disallowFalseMove)
                {
                    if (res && !unsetAction.Successful)
                        actionsPerformed.Add(unsetAction);
                    Unperform();
                    return false;
                }
                else if (!res || !unsetAction.Successful)
                {
                    successful = false;
                }
                actionsPerformed.Add(unsetAction);
            }

            if (targetState != EdgeState.Empty)
            {
                bool res = mesh.Perform(edgeIndex, targetState, actionsPerformed, 0);
                if (!res && disallowFalseMove)
                {
                    successful = false;
                    Unperform();
                    return false;
                }
                else if (!res)
                {
                    successful = false;
                }
            }

            return true;
        }

        public void Unperform()
        {
            if (actionsPerformed != null && actionsPerformed.Count > 0)
            {
                mesh.Unperform(actionsPerformed);
            }
        }

        public bool Equals(IAction? other)
        {
            return other is PuzzleSetEdgeStateAction o && o.mesh == mesh && o.edgeIndex == edgeIndex && o.targetState == targetState && o.disallowFalseMove == disallowFalseMove;
        }
    }
}
