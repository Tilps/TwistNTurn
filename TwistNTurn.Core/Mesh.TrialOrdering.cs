using System.Collections.Generic;

namespace TwistNTurn
{
    public enum PostProgressAction
    {
        Continue,
        ReorderRemaining,
        RestartPass,
    }

    public interface ITrialOrderingPolicy
    {
        /// <summary>
        /// Called once at the beginning of TrySolve after PerformStart completes.
        /// Policies with cached per-solve state should reset caches here.
        /// </summary>
        void OnSolveStart(Mesh mesh) { }

        /// <summary>
        /// Called at the start of IterativeTrySolveInternal after gathering trials.
        /// May reorder trials[startIndex..] in place.
        /// </summary>
        void OrderIterativeTrials(Mesh mesh, List<IAction> trials, int startIndex, bool isTopLevel, int iterativeRecDepth);

        /// <summary>
        /// Called after a trial commits real progress in IterativeTrySolveInternal.
        /// </summary>
        PostProgressAction OnTrialProgress(Mesh mesh, List<IAction> trials, int currentIndex, int addedActions, bool isTopLevel, int iterativeRecDepth)
            => PostProgressAction.Continue;

        /// <summary>
        /// Optional: choose index in [startIndex, trials.Count) for RecursiveTrySolveInternal to branch on next.
        /// Return -1 to preserve default sequential order.
        /// </summary>
        int SelectRecursiveTrial(Mesh mesh, List<IAction> trials, int startIndex)
            => -1;
    }

    public partial class Mesh
    {
        /// <summary>
        /// Optional trial ordering policy used by IterativeTrySolve / RecursiveTrySolve.
        /// Null preserves exact upstream sequential behavior with zero overhead.
        /// Note: stateful policies (such as IncidenceGnnPolicy) hold per-solve scratch buffers
        /// and should not be shared concurrently across threads.
        /// </summary>
        public ITrialOrderingPolicy? TrialOrderingPolicy { get; set; }

        public IReadOnlyList<List<int>> EdgeSets => edgeSets;
        public IReadOnlyList<List<int>> ColorSets => colorSets;
        public IReadOnlyList<List<int>> CellColorSets => cellColorSets;

        public bool IsPointlessTrialAction(IAction action) => IsPointlessTrial(action);

        public EdgePairRestriction GetEdgePairRestriction(int e1, int e2)
        {
            return edgePairRestrictions != null ? edgePairRestrictions[e1, e2] : EdgePairRestriction.None;
        }

        public int CountEmptyEdges()
        {
            int count = 0;
            for (int i = 0; i < edges.Count; i++)
            {
                if (edges[i].State == EdgeState.Empty)
                    count++;
            }
            return count;
        }
    }
}
