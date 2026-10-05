using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace TwistNTurn
{
    public partial class Mesh
    {
        public int IterativeSolverDepth
        {
            get
            {
                return iterativeSolverDepth;
            }
            set
            {
                iterativeSolverDepth = value;
            }
        }
        private int iterativeSolverDepth = 1;

        public int IterativeRecMaxDepth
        {
            get
            {
                return iterativeRecMaxDepth;
            }
            set
            {
                iterativeRecMaxDepth = value;
            }
        }
        private int iterativeRecMaxDepth = 2;
        private int iterativeRecDepth = 1;


        public bool UseMerging
        {
            get
            {
                return useMerging;
            }
            set
            {
                useMerging = value;
            }
        }
        private bool useMerging;

        public bool ConsiderMultipleLoops
        {
            get
            {
                return considerMultipleLoops;
            }
            set
            {
                considerMultipleLoops = value;
            }
        }
        private bool considerMultipleLoops = false;

        public Mesh? FinalSolution
        {
            get
            {
                return finalSolution;
            }
        }
        Mesh? finalSolution;

        public int[]? FinalDepthPatern
        {
            get
            {
                return finalDepthPatern;
            }
        }
        int[]? finalDepthPatern;

        public Mesh SolutionFound
        {
            get
            {
                return solutionsFound[0];
            }
        }
        List<Mesh> solutionsFound = new List<Mesh>();

        public int[] DepthPatern
        {
            get
            {
                return solutionDepthPatern[0];
            }
        }
        List<int[]> solutionDepthPatern = new List<int[]>();
        List<int> curDepthPatern = new List<int>();


        public SolveState TrySolve()
        {
            return TrySolve(false);
        }

        public SolveState TrySolve(bool noRollback)
        {
            solutionsFound.Clear();
            solutionDepthPatern.Clear();
            List<IAction> edgeChanges = new List<IAction>();
            try
            {
                if (!PerformStart(edgeChanges))
                    return SolveState.NoSolutions;
                // No point doing anything if we're already complete.
                PerformEndIfPossible(edgeChanges);
                curDepthPatern.Clear();
                if (solverMethod == SolverMethod.Iterative)
                    return IterativeTrySolve(noRollback);
                else
                    return RecursiveTrySolve(noRollback);
            }
            catch (OperationCanceledException)
            {
                // We're screwed, no point trying to rollback changes.
                edgeChanges.Clear();
                return SolveState.NoSolutions;
            }
            finally
            {
                if (!noRollback)
                    Unperform(edgeChanges);
            }
        }

        private bool topLevel = false;

        public bool PerformStart(List<IAction> changes)
        {
            try
            {
                topLevel = true;
                if (superSlowMo)
                {
                    for (int i = 0; i < intersections.Count; i++)
                        if (intersections[i].Type != IntersType.Unknown)
                            if (!Perform(null, changes, int.MaxValue, new List<int> { i }))
                                return false;
                    return true;
                }
                else
                {
                    List<int> allIntersWithType = new List<int>();
                    for (int i = 0; i < intersections.Count; i++)
                        if (intersections[i].Type != IntersType.Unknown)
                            allIntersWithType.Add(i);
                    // TODO: turn on coloring cheats around this?
                    return Perform(null, changes, int.MaxValue, allIntersWithType);
                }
            }
            finally
            {
                topLevel = false;
            }
        }

        /// <summary>
        /// Use iterative solver on full power before using recursion to finish off if needed.
        /// </summary>
        public bool ContaminateFullSolver
        {
            get
            {
                return contaminateFullSolver;
            }
            set
            {
                contaminateFullSolver = value;
            }
        }
        private bool contaminateFullSolver = true;

        private SolveState RecursiveTrySolve(bool noRollback)
        {
            List<IAction> realChanges = new List<IAction>();
            int backupDepth = iterativeSolverDepth;
            try
            {
                if (contaminateFullSolver)
                {
                    iterativeSolverDepth = int.MaxValue;
                    SolveState iterativeTry = IterativeTrySolveWithoutRollback(realChanges);
                    if (iterativeTry == SolveState.Solved)
                        return SolveState.Solved;
                    else if (iterativeTry == SolveState.NoSolutions)
                        return SolveState.NoSolutions;
                    iterativeSolverDepth = backupDepth;
                }
                if (MeshChangeUpdate != null)
                    MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, null, true, true));
                List<IAction> trials = new List<IAction>();
                if (useCellColoring && useCellColoringTrials)
                    GatherAllCellTrials(trials);
                else
                    GatherAll(trials);
                return RecursiveTrySolveInternal(trials, 0);
            }
            catch (OperationCanceledException)
            {
                // We're screwed, no point trying to rollback changes.
                realChanges.Clear();
                return SolveState.NoSolutions;
            }
            finally
            {
                iterativeSolverDepth = backupDepth;
                if (!noRollback)
                    Unperform(realChanges);
            }
        }

        private SolveState RecursiveTrySolveInternal(List<IAction> trials, int index)
        {
            if (pruning && earlyFail)
            {
                earlyFail = false;
                return SolveState.MultipleSolutions;
            }
            // Multiple solutions, we can make progress from here using recursion.
            List<IAction> edgeChanges1 = new List<IAction>();
            List<IAction> edgeChanges2 = new List<IAction>();
            while (index < trials.Count && IsPointlessTrial(trials[index]))
                index++;
            curDepthPatern.Add(index);
            if (index < trials.Count)
            {
                bool success1 = Perform(trials[index], edgeChanges1);
                if (MeshChangeUpdate != null)
                    MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, null, true, true));
                if (success1)
                {
                    SolveState result = RecursiveTrySolveInternal(trials, index + 1);
                    if (result == SolveState.MultipleSolutions)
                    {
                        Unperform(edgeChanges1);
                        curDepthPatern.RemoveAt(curDepthPatern.Count - 1);
                        return SolveState.MultipleSolutions;
                    }
                    else if (result == SolveState.NoSolutions)
                        success1 = false;
                }
                Unperform(edgeChanges1);
                bool success2 = Perform(GetOppositeAction(trials[index]), edgeChanges2);
                if (MeshChangeUpdate != null)
                    MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, null, true, true));
                if (success2)
                {
                    SolveState result = RecursiveTrySolveInternal(trials, index + 1);
                    if (result == SolveState.MultipleSolutions)
                    {
                        Unperform(edgeChanges2);
                        curDepthPatern.RemoveAt(curDepthPatern.Count - 1);
                        return SolveState.MultipleSolutions;
                    }
                    else if (result == SolveState.NoSolutions)
                        success2 = false;
                }
                Unperform(edgeChanges2);
                curDepthPatern.RemoveAt(curDepthPatern.Count - 1);
                if (!success1 && !success2)
                    return SolveState.NoSolutions;
                else if (success1 && success2)
                    return SolveState.MultipleSolutions;
                else
                    return SolveState.Solved;
            }
            solutionsFound.Add(this.Clone());
            solutionDepthPatern.Add(curDepthPatern.ToArray());
            curDepthPatern.RemoveAt(curDepthPatern.Count - 1);
            if (solutionsFound.Count > 1)
                return SolveState.MultipleSolutions;
            return SolveState.Solved;
        }

        public double PercentSolved
        {
            get
            {
                return percentSolved;
            }
        }
        private double percentSolved = 0.0;

        private SolveState IterativeTrySolve(bool noRollback)
        {
            percentSolved = 0.0;
            List<IAction> realChanges = new List<IAction>();
            try
            {
                SolveState res = IterativeTrySolveWithoutRollback(realChanges);
                int count = 0;
                foreach (Edge edge in edges)
                {
                    if (edge.State == EdgeState.Empty)
                        count++;
                }
                percentSolved = (double)count / (double)edges.Count;
                return res;
            }
            catch (OperationCanceledException)
            {
                // We're screwed anyway, no point trying to undo the changes.
                realChanges.Clear();
                return SolveState.NoSolutions;
            }
            finally
            {
                if (!noRollback)
                    Unperform(realChanges);
            }
        }

        private SolveState IterativeTrySolveWithoutRollback(List<IAction> realChanges)
        {
            int count = 0;
            bool changed = true;
            iterativeRecDepth = 1;
            if (pruning)
            {
                oldFirstOrder = null;
            }
            while (changed)
            {
                count = 0;
                foreach (Edge edge in edges)
                {
                    if (edge.State == EdgeState.Empty)
                        count++;
                }
                curDepthPatern.Add(count);
                changed = false;
                if (!IterativeTrySolveInternal(realChanges, ref changed, iterativeRecDepth, null))
                {
                    return SolveState.NoSolutions;
                }
                if (changed == false)
                {
                    changed = PerformEndIfPossible(realChanges);
                }
                if (!changed)
                {
                    if (iterativeRecDepth < iterativeRecMaxDepth)
                    {
                        iterativeRecDepth++;
                        changed = true;
                    }
                }
                if (changed && !Partition())
                    return SolveState.MultipleSolutions;
            }
            if (count > 0)
                return SolveState.MultipleSolutions;
            solutionsFound.Add(this.Clone());
            solutionDepthPatern.Add(curDepthPatern.ToArray());
            return SolveState.Solved;
        }

        private bool Partition()
        {/* Failed to provide noticible performance improvement, not worth doing until we use it to decompose the puzzle and actually help solve things rather than just to detect non-progress faster.
            if (pruning)
            {
                // Retrieve existing partition sizes.
                // Calculate new partitions.
                // If any first order partition has not changed size, return false.
                // We can't do any second order partitions unless the solver understands the implicit edgeset link from the first order partitions.
                DisjointTracker tracker = new DisjointTracker(edges.Count);
                for (int i = 0; i < edges.Count; i++)
                {
                    if (edges[i].State == EdgeState.Empty)
                    {
                        tracker.Add(i);
                    }
                }
                for (int i = 0; i < edges.Count; i++)
                {
                    if (edges[i].State == EdgeState.Empty)
                    {
                        List<int> affecting = GetAffectingEdges(i);
                        foreach (int a in affecting)
                        {
                            if (edges[a].State == EdgeState.Empty)
                            {
                                tracker.Union(i, a);
                            }
                        }
                    }
                }
                Dictionary<int, List<int>> meetPoints = new Dictionary<int, List<int>>();
                List<int> representatives = new List<int>();
                for (int i = 0; i < intersections.Count; i++)
                {
                    // at each intersection check edges to see if they come from differing sets, if so then it is a meet point.
                    // If a set has more than 2 meetPoints it is not first order.  If it has 0 meet points we additionally do not care.
                    // Need to follow edge sets as well.
                    List<int> realIs = new List<int>();
                    realIs.Add(i);
                    Intersection inters = intersections[i];
                    if (inters.EdgeSet >= 0)
                    {
                        if (inters.FilledCount == 1)
                        {
                            // Follow filled lines to get other end.
                            int lastEdge = GetNextEdge(i, -1);
                            int nextI = GetOtherInters(lastEdge, i);
                            while (intersections[nextI].FilledCount > 1)
                            {
                                int nextEdge = GetNextEdge(nextI, lastEdge);
                                nextI = GetOtherInters(nextEdge, nextI);
                                lastEdge = nextEdge;
                            }
                            if (nextI < i)
                                continue;
                            realIs.Add(nextI);
                        }

                    }
                    representatives.Clear();
                    foreach (int realI in realIs)
                    {
                        foreach (int e1 in intersections[realI].Edges)
                        {
                            if (edges[e1].State == EdgeState.Empty)
                            {
                                int rep = tracker.GetRepresentative(e1);
                                if (!representatives.Contains(rep))
                                    representatives.Add(rep);
                            }
                        }
                    }
                    if (representatives.Count > 1)
                    {
                        for (int j = 0; j < representatives.Count; j++)
                        {
                            List<int> points;
                            if (!meetPoints.TryGetValue(representatives[j], out points))
                            {
                                points = new List<int>();
                                meetPoints.Add(representatives[j], points);
                            }
                            points.Add(i);
                        }
                    }
                }
                Dictionary<int, int> firstOrder = new Dictionary<int, int>();
                foreach (KeyValuePair<int, List<int>> kvp in meetPoints)
                {
                    if (kvp.Value.Count == 2)
                        firstOrder[kvp.Key] = 0;
                }
                for (int i = 0; i < edges.Count; i++)
                {
                    if (edges[i].State == EdgeState.Empty)
                    {
                        int rep = tracker.GetRepresentative(i);
                        int value;
                        if (firstOrder.TryGetValue(rep, out value))
                        {
                            firstOrder[rep] = value + 1;
                        }
                    }
                }
                if (oldFirstOrder != null)
                {
                    foreach (KeyValuePair<int, int> kvp in oldFirstOrder)
                    {
                        if (edges[kvp.Key].State == EdgeState.Empty)
                        {
                            int rep = tracker.GetRepresentative(kvp.Key);
                            int value;
                            if (firstOrder.TryGetValue(rep, out value))
                            {
                                if (value == kvp.Value)
                                    return false;
                            }
                        }
                    }
                }
                oldFirstOrder = firstOrder;
            }
          */
            return true;
        }
        private Dictionary<int, int>? oldFirstOrder;
        /*
        private List<int> GetAffectingEdges(int i)
        {
            // TODO: this approach fails to recognize potential second order effects such as 'the airlock'
            // consider an area which touches outside areas at multiple locations, with locked pairs.  It can be trivially seen that if there are numbers 
            // outside of the area, the entire area must be empty.  Thus the outside areas may be independent, but will seem connected via the airlock.
            // We may need a seperate pass to detect and fill those so as to allow this to be effective.
            List<int> result = new List<int>();
            foreach (int cellIndex in edges[i].Cells)
            {
                if (cells[cellIndex].TargetCount >= 0)
                {
                    // TODO: optimize based on antilock edges causing seperation.
                    foreach (int e1 in cells[cellIndex].Edges)
                    {
                        if (e1 != i && !result.Contains(e1))
                            result.Add(e1);
                    }
                }
            }
            foreach (int intersIndex in edges[i].Intersections)
            {
                // TODO: optimize based on antilock edges causing some edges to be irrelivent.
                foreach (int e1 in intersections[intersIndex].Edges)
                {
                    if (e1 != i && !result.Contains(e1))
                        result.Add(e1);
                }
            }
            return result;
        }*/

#if DEBUG
        private void ValidateMaximalProgression()
        {
            // TODO: add more maximal progression tests.
            if (useCellColoring)
            {
                foreach (Cell c in cells)
                {
                    if (Math.Abs(c.Color) != 1)
                        continue;
                    foreach (int edgeIndex in c.Edges)
                    {
                        Edge edge = edges[edgeIndex];
                        bool found = false;
                        foreach (int cellIndex in edge.Cells)
                        {
                            if (cells[cellIndex] == c)
                                continue;
                            found = true;
                            if (Math.Abs(cells[cellIndex].Color) == 1)
                            {
                                if (edge.State == EdgeState.Empty)
                                    throw new Exception("Solver failed basic consistancy, empty edge.");
                            }
                            else
                            {
                                if (edge.State != EdgeState.Empty)
                                    throw new Exception("Solver failed basic consistancy, cell color missing.");
                            }
                        }
                        if (!found)
                        {
                            if (edge.State == EdgeState.Empty)
                                throw new Exception("Solver failed basic consistancy, empty outside edge.");
                        }
                    }
                }
            }
        }
#endif

        private bool IterativeTrySolveInternal(List<IAction> realChanges, ref bool changed, int iterativeRecDepth, IAction? locusAction)
        {
            if (iterativeSolverDepth < 0)
                return true;
            bool iterTopLevel = true;
            List<IAction> trials = new List<IAction>();
            if (iterativeRecDepth != this.iterativeRecDepth)
            {
                iterTopLevel = false;
                if (locusAction != null) GatherLocals(realChanges, trials, locusAction);
            }
            else
            {
                GatherAll(trials);
            }
            List<IAction> edgeChanges1 = new List<IAction>();
            List<IAction> edgeChanges2 = new List<IAction>();
            IAction? lastTrial = null;
            for (int i = 0; i < trials.Count; i++)
            {
                if (pruning && earlyFail)
                {
                    if (iterTopLevel)
                    {
                        earlyFail = false;
                    }
                    return false; // This will incorrectly mark as no-solutions rather than multiple solutions... but generator doesn't care at the moment.
                    // TODO: provide mechanism to early out with multiple solutions from the iterative solver.
                }
                if (lastTrial != null && lastTrial.Equals(trials[i]))
                    continue;
                lastTrial = trials[i];
                if (IsPointlessTrial(trials[i]))
                    continue;
                edgeChanges1.Clear();
                edgeChanges2.Clear();
                bool success1 = Perform(trials[i], edgeChanges1, iterativeSolverDepth);
                if (success1 && iterativeRecDepth > 1)
                {
                    bool ignored = false;
                    success1 = IterativeTrySolveInternal(edgeChanges1, ref ignored, iterativeRecDepth - 1, trials[i]);
                }
                if (MeshChangeUpdate != null)
                    MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, edgeChanges1, success1));
                Unperform(edgeChanges1);
                bool success2 = Perform(GetOppositeAction(trials[i]), edgeChanges2, iterativeSolverDepth);
                if (success2 && iterativeRecDepth > 1)
                {
                    bool ignored = false;
                    success2 = IterativeTrySolveInternal(edgeChanges2, ref ignored, iterativeRecDepth - 1, trials[i]);
                }
                if (MeshChangeUpdate != null)
                    MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, edgeChanges2, success2));
                Unperform(edgeChanges2);
                if (!success1 && !success2)
                    return false;
                else if (success1 && success2)
                {
                    List<IAction> merged = new List<IAction>();
                    if (UseColoring && UseDerivedColoring)
                    {
                        merged.AddRange(DeriveColoring(edgeChanges1, edgeChanges2));
                    }
                    // TODO: consider derived edge restrictions, although I think merging is probably sufficient.
                    if (useMerging)
                    {
                        merged.AddRange(MergeChanges(edgeChanges1, edgeChanges2));
                    }
                    if (merged.Count > 0)
                    {
                        bool colorCheatBackup = coloringCheats;
                        try
                        {
                            if (iterTopLevel)
                            {
                                coloringCheats = true;
                                topLevel = true;
                            }
                            int length = realChanges.Count;
                            if (!Perform(merged, realChanges, iterTopLevel))
                                return false;

                            if (MeshChangeUpdate != null)
                                MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, null, true));
                            if (realChanges.Count != length)
                                changed = true;
                            // This doesn't appear to make any sense and PerformStart has changed now.
                            /*if (changed && topLevel)
                            {
                                if (!PerformStart(realChanges))
                                    return false;
                            }*/
                        }
                        finally
                        {
                            coloringCheats = colorCheatBackup;
                            topLevel = false;
                        }
                    }
                }
                else if (success1)
                {
                    bool colorCheatBackup = coloringCheats;
                    try
                    {
                        if (iterTopLevel)
                        {
                            coloringCheats = true;
                            topLevel = true;
                        }
                        int length = realChanges.Count;
                        if (!Perform(edgeChanges1, realChanges, iterTopLevel))
                            return false;
                        if (MeshChangeUpdate != null)
                            MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, null, true));
                        if (realChanges.Count != length)
                            changed = true;
                    }
                    finally
                    {
                        coloringCheats = colorCheatBackup;
                        topLevel = false;
                    }
                }
                else
                {
                    bool colorCheatBackup = coloringCheats;
                    try
                    {
                        if (iterTopLevel)
                        {
                            coloringCheats = true;
                            topLevel = true;
                        }
                        int length = realChanges.Count;
                        if (!Perform(edgeChanges2, realChanges, iterTopLevel))
                            return false;
                        if (MeshChangeUpdate != null)
                            MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, null, true));
                        if (realChanges.Count != length)
                            changed = true;
                    }
                    finally
                    {
                        coloringCheats = colorCheatBackup;
                        topLevel = false;
                    }
                }
            }
            return true;
        }

        private void GatherAll(List<IAction> trials)
        {
            for (int i = 0; i < edges.Count; i++)
            {
                trials.Add(new SetAction(this, i, EdgeState.Filled));
            }
            GatherAllCellTrials(trials);
        }

        private void GatherAllCellTrials(List<IAction> trials)
        {
            if (useCellColoring && useCellColoringTrials)
            {
                for (int i = 0; i < cells.Count; i++)
                {
                    trials.Add(new CellColorJoinAction(this, i, -1, true));
                }
            }
        }

        private IAction GetOppositeAction(IAction iAction)
        {
            if (iAction is SetAction)
            {
                SetAction setAction = (SetAction)iAction;
                return new SetAction(this, setAction.EdgeIndex, FlipEdgeState(setAction.EdgeState));
            }
            else if (iAction is ColorJoinAction)
            {
                ColorJoinAction cjAction = (ColorJoinAction)iAction;
                return new ColorJoinAction(this, cjAction.Edge1, cjAction.Edge2, !cjAction.Same);
            }
            else if (iAction is CellColorJoinAction)
            {
                CellColorJoinAction ccjAction = (CellColorJoinAction)iAction;
                return new CellColorJoinAction(this, ccjAction.Cell1, ccjAction.Cell2, !ccjAction.Same);
            }
            throw new NotSupportedException("GetOppositeAction doesn't support the passed action type.");
        }

        private EdgeState FlipEdgeState(EdgeState edgeState)
        {
            if (edgeState == EdgeState.Filled)
                return EdgeState.Excluded;
            if (edgeState == EdgeState.Excluded)
                return EdgeState.Filled;
            throw new InvalidOperationException("Cannot flip an empty edge.");
        }

        private bool IsPointlessTrial(IAction iAction)
        {
            if (iAction is SetAction)
            {
                SetAction setAction = (SetAction)iAction;
                return edges[setAction.EdgeIndex].State != EdgeState.Empty;
            }
            else if (iAction is ColorJoinAction)
            {
                ColorJoinAction cjAction = (ColorJoinAction)iAction;
                return Math.Abs(edges[cjAction.Edge1].Color) == Math.Abs(edges[cjAction.Edge2].Color);
            }
            else if (iAction is CellColorJoinAction)
            {
                CellColorJoinAction ccjAction = (CellColorJoinAction)iAction;
                if (ccjAction.Cell2 == -1)
                    return Math.Abs(cells[ccjAction.Cell1].Color) == 1;
                else
                    return Math.Abs(cells[ccjAction.Cell1].Color) == Math.Abs(cells[ccjAction.Cell2].Color);
            }
            throw new NotSupportedException("IsPointless doesn't support the passed action type.");

        }

        private IEnumerable<IAction> DeriveColoring(List<IAction> edgeChanges1, List<IAction> edgeChanges2)
        {
            List<IAction> derived = new List<IAction>();
            EdgeState[] side1 = new EdgeState[edges.Count];
            EdgeState[] side2 = new EdgeState[edges.Count];
            foreach (IAction change in edgeChanges1)
            {
                if (change is SetAction)
                {
                    SetAction sa = (SetAction)change;
                    side1[sa.EdgeIndex] = sa.EdgeState;
                }
            }
            foreach (IAction change in edgeChanges2)
            {
                if (change is SetAction)
                {
                    SetAction sa = (SetAction)change;
                    side2[sa.EdgeIndex] = sa.EdgeState;
                }
            }
            int firstI = -1;
            for (int i = 0; i < edges.Count; i++)
            {
                if (side2[i] != side1[i] && side2[i] != EdgeState.Empty && side1[i] != EdgeState.Empty)
                {
                    if (firstI == -1)
                        firstI = i;
                    else
                        derived.Add(new ColorJoinAction(this, firstI, i, side1[i] == side1[firstI]));
                }
            }
            return derived;
        }

        private void GatherLocals(List<IAction> realChanges, List<IAction> trials, IAction locusAction)
        {
            for (var i = 0; i < realChanges.Count; i++)
            {
                IAction action = realChanges[i];
                if (action is SetAction)
                {
                    SetAction setAction = (SetAction) action;
                    for (var j = 0; j < setAction.GetAffectedEdges().Count; j++)
                    {
                        int edge = setAction.GetAffectedEdges()[j];
                        GatherNearEdge(trials, edge, locusAction, true);
                    }
                }
                else if (action is ColorJoinAction)
                {
                    ColorJoinAction cjAction = (ColorJoinAction) action;
                    for (var j = 0; j < cjAction.GetAffectedEdges().Count; j++)
                    {
                        int edge = cjAction.GetAffectedEdges()[j];
                        GatherNearEdge(trials, edge, locusAction, true);
                    }
                }
                else if (action is CellColorJoinAction)
                {
                    CellColorJoinAction ccjAction = (CellColorJoinAction) action;
                    for (var j = 0; j < ccjAction.GetAffectedCells().Count; j++)
                    {
                        int cell = ccjAction.GetAffectedCells()[j];
                        for (var k = 0; k < cells[cell].Edges.Count; k++)
                        {
                            int edge = cells[cell].Edges[k];
                            GatherNearEdge(trials, edge, locusAction, false);
                        }
                    }
                }
            }
            trials.Sort(new ActionSorter());
        }

        private void GatherNearEdge(List<IAction> smarts, int edge, IAction locusAction, bool followCells)
        {
            int maxDist = iterativeSolverDepth > edges.Count ? iterativeSolverDepth : (iterativeSolverDepth + 2); // TODO: this should really be + half the max edge count for any cell.

            // This gets the same edge upto 4 times, so uniquification is probably useful by the caller (faster then us doing it.)
            Edge e = edges[edge];
            for (var index = 0; index < intersections[e.Intersections[0]].Edges.Count; index++)
            {
                int i = intersections[e.Intersections[0]].Edges[index];
                if (edges[i].State == EdgeState.Empty && GetEdgeDistance(i, locusAction) <= maxDist)
                    smarts.Add(new SetAction(this, i, EdgeState.Filled));
            }
            for (var index = 0; index < intersections[e.Intersections[1]].Edges.Count; index++)
            {
                int i = intersections[e.Intersections[1]].Edges[index];
                if (edges[i].State == EdgeState.Empty && GetEdgeDistance(i, locusAction) <= maxDist)
                    smarts.Add(new SetAction(this, i, EdgeState.Filled));
            }
            if (useCellColoring && useCellColoringTrials)
            {
                for (var index = 0; index < e.Cells.Count; index++)
                {
                    int i = e.Cells[index];
                    if (Math.Abs(cells[i].Color) != 1 && GetCellDistance(i, locusAction) <= maxDist)
                        smarts.Add(new CellColorJoinAction(this, i, -1, true));
                }
            }
            if (followCells)
            {
                for (var index = 0; index < cells[e.Cells[0]].Edges.Count; index++)
                {
                    int i = cells[e.Cells[0]].Edges[index];
                    if (edges[i].State == EdgeState.Empty && GetEdgeDistance(i, locusAction) <= maxDist)
                        smarts.Add(new SetAction(this, i, EdgeState.Filled));
                }
                if (e.Cells.Count > 1)
                {
                    for (var index = 0; index < cells[e.Cells[1]].Edges.Count; index++)
                    {
                        int i = cells[e.Cells[1]].Edges[index];
                        if (edges[i].State == EdgeState.Empty && GetEdgeDistance(i, locusAction) <= maxDist)
                            smarts.Add(new SetAction(this, i, EdgeState.Filled));
                    }
                }
            }
        }

        private int GetCellDistance(int cell, IAction locusAction)
        {
            List<int> locusEdges = GetLocusEdges(locusAction);
            List<int> cellEdges = cells[cell].Edges;
            int minDist = int.MaxValue;
            for (int i = 0; i < locusEdges.Count; i++)
            {
                for (int j = 0; j < cellEdges.Count; j++)
                {
                    int dist = edgeDistances[locusEdges[i], cellEdges[j]];
                    if (dist < minDist)
                        minDist = dist;
                }
            }
            return minDist;
        }

        private List<int> GetLocusEdges(IAction locusAction)
        {
            List<int> result = new List<int>();
            if (locusAction is SetAction)
            {
                SetAction setAction = (SetAction)locusAction;
                result.Add(setAction.EdgeIndex);
                return result;
            }
            else if (locusAction is ColorJoinAction)
            {
                ColorJoinAction cjAction = (ColorJoinAction)locusAction;
                result.Add(cjAction.Edge1);
                result.Add(cjAction.Edge2);
                return result;
            }
            else if (locusAction is CellColorJoinAction)
            {
                CellColorJoinAction ccjAction = (CellColorJoinAction)locusAction;
                result.AddRange(cells[ccjAction.Cell1].Edges);
                if (ccjAction.Cell2 != -1)
                    result.AddRange(cells[ccjAction.Cell2].Edges);
                return result;
            }
            throw new NotSupportedException("GetLocusEdges doesn't support the passed action type.");
        }

        private int GetEdgeDistance(int edge, IAction locusAction)
        {
            List<int> locusEdges = GetLocusEdges(locusAction);
            int minDist = int.MaxValue;
            for (int i = 0; i < locusEdges.Count; i++)
            {
                int dist = edgeDistances[locusEdges[i], edge];
                if (dist < minDist)
                    minDist = dist;
            }
            return minDist;
        }

        private bool PerformEndIfPossible(List<IAction> realChanges)
        {
            bool changed = false;
            if (satisifiedIntersCount == intersections.Count)
            {
                int edgeSetCount = 0;
                foreach (List<int> edgeSet in edgeSets)
                {
                    if (edgeSet.Count > 0)
                        edgeSetCount++;
                }
                // If there is exactly 1 edge set and all intersections are happy, we are in a complete situation.
                if (edgeSetCount == 1)
                {
                    for (int i = 0; i < edges.Count; i++)
                    {
                        Edge e = edges[i];
                        if (e.State == EdgeState.Empty)
                        {
                            changed = true;
                            Perform(i, EdgeState.Excluded, realChanges, 0);
                        }
                    }
                }
            }
            return changed;
        }

        private Mesh Clone()
        {
            return new Mesh(this);
        }
        public bool PerformListNoRecurse(List<IAction> list)
        {
            foreach (IAction action in list)
            {
                if (!action.Perform() || !action.Successful)
                    return false;
            }
            return true;

        }
        public void PerformListRegardless(List<IAction> list)
        {
            foreach (IAction action in list)
            {
                action.Perform();
            }

        }
        private bool Perform(List<IAction> list, List<IAction> realChanges, bool useFullPower)
        {
            foreach (IAction action in list)
            {
                if (action is SetAction)
                {
                    SetAction act = (SetAction)action;
                    Edge edge = edges[act.EdgeIndex];
                    if (edge.State == EdgeState.Empty)
                    {
                        if (!Perform(act.EdgeIndex, act.EdgeState, realChanges, useFullPower ? int.MaxValue : iterativeSolverDepth))
                            return false;
                    }
                    else if (edge.State != act.EdgeState)
                        throw new InvalidOperationException("doh");
                }
                else if (action is ColorJoinAction)
                {
                    if (!action.Perform())
                        return false;
                    if (!action.Successful)
                        return false;
                    if (!((ColorJoinAction)action).WasteOfTime)
                        realChanges.Add(action);
                }
                else if (action is CellColorJoinAction)
                {
                    if (!action.Perform())
                        return false;
                    if (!action.Successful)
                        return false;
                    if (!((CellColorJoinAction)action).WasteOfTime)
                        realChanges.Add(action);
                }
                else if (action is EdgeRestrictionAction)
                {
                    if (!action.Perform())
                        return false;
                    if (!action.Successful)
                        return false;
                    if (!((EdgeRestrictionAction)action).WasteOfTime)
                        realChanges.Add(action);
                }
                else
                    throw new NotSupportedException("Perform List only supports SetAction and colorjoin action.");
            }
            return true;
        }

        private List<IAction> MergeChanges(List<IAction> edgeChanges1, List<IAction> edgeChanges2)
        {
            List<IAction> res = new List<IAction>();
            if (edgeChanges2.Count * edgeChanges1.Count > 1000)
            {
                Dictionary<IAction, bool> lookup = new Dictionary<IAction, bool>();
                foreach (IAction action in edgeChanges2)
                    lookup[action] = true;
                foreach (IAction action in edgeChanges1)
                    if (lookup.ContainsKey(action))
                        res.Add(action);
            }
            else
            {
                foreach (IAction action in edgeChanges1)
                    if (edgeChanges2.Contains(action))
                        res.Add(action);
            }
            return res;
        }

        public void Unperform(List<IAction> backup)
        {
            for (int i = backup.Count - 1; i >= 0; i--)
            {
                IAction backupBit = backup[i];
                Unperform(backupBit);
            }
        }

        public void Unperform(IAction backup)
        {
            backup.Unperform();
        }

        public bool Perform(int edgeIndex, EdgeState state, List<IAction> backup)
        {
            return Perform(edgeIndex, state, backup, int.MaxValue);
        }

        EdgeState[] edgesSeen = null!;
        TriState[,] edgePairsSeen = null!;
        List<KeyValuePair<int, int>> edgePairsToClean = new List<KeyValuePair<int, int>>();
        TriState[,] cellPairsSeen = null!;
        List<KeyValuePair<int, int>> cellPairsToClean = new List<KeyValuePair<int, int>>();
        EdgePairRestriction[,] edgeRestrictsSeen = null!;
        List<KeyValuePair<int, int>> edgeRestrictsToClean = new List<KeyValuePair<int, int>>();
        bool[] cellsSeen = null!;
        bool[] cellColorEdgeColorsSeen = null!;
        bool[] intersectsSeen = null!;
        bool[] interactsSeen = null!;
        bool[] interactsSeen2 = null!;
        bool[] interactsSeen3 = null!;
        private bool[] colorSetsSeen = null!;

        public bool UseColoring
        {
            get
            {
                return useColoring;
            }
            set
            {
                useColoring = value;
            }
        }
        private bool useColoring = false;

        public bool UseCellPairs
        {
            get
            {
                return useCellPairs;
            }
            set
            {
                useCellPairs = value;
            }
        }
        private bool useCellPairs = false;

        public bool UseCellPairsTopLevel
        {
            get
            {
                return useCellPairsTopLevel;
            }
            set
            {
                useCellPairsTopLevel = value;
            }
        }
        private bool useCellPairsTopLevel = false;

        public bool UseEdgeRestricts
        {
            get
            {
                return useEdgeRestricts;
            }
            set
            {
                useEdgeRestricts = value;
            }
        }
        private bool useEdgeRestricts = false;

        public bool UseDerivedColoring
        {
            get
            {
                return useDerivedColoring;
            }
            set
            {
                useDerivedColoring = value;
            }
        }
        private bool useDerivedColoring = true;

        public bool UseCellColoring
        {
            get
            {
                return useCellColoring;
            }
            set
            {
                useCellColoring = value;
            }
        }
        private bool useCellColoring = false;

        public bool UseCellColoringTrials
        {
            get
            {
                return useCellColoringTrials;
            }
            set
            {
                useCellColoringTrials = value;
            }
        }
        private bool useCellColoringTrials = true;

        public bool ColoringCheats
        {
            get
            {
                return coloringCheats;
            }
            set
            {
                coloringCheats = value;
            }
        }
        private bool coloringCheats = false;

        public bool SuperSlowMo
        {
            get
            {
                return superSlowMo;
            }
            set
            {
                superSlowMo = value;
            }
        }
        private bool superSlowMo = false;

        public bool Perform(int edgeIndex, EdgeState state, List<IAction> backup, int maxDepth)
        {
            return Perform(new SetAction(this, edgeIndex, state), backup, maxDepth, null);
        }

        public int deepestNonEmpty = -1;
        private bool Perform(IAction? action, List<IAction> backup)
        {
            return Perform(action, backup, int.MaxValue);
        }
        private bool Perform(IAction? action, List<IAction> backup, int maxDepth)
        {
            return Perform(action, backup, maxDepth, null);
        }

        List<IAction>[] moves = null!;
        List<int>[] toConsiderEdges = null!;
        List<int>[] toConsiderEdgeSets = null!;
        List<int>[] toConsiderEdgeColors = null!;
        List<int>[] toConsiderIntersects = null!;
        List<int>[] toConsiderCellColors = null!;

        private void AllocateForPerform(int size)
        {
            if (moves == null || size > moves.Length)
            {
                moves = new List<IAction>[size];
                for (int i = 0; i < moves.Length; i++)
                    moves[i] = new List<IAction>();

                toConsiderEdges = new List<int>[size];
                for (int i = 0; i < size; i++)
                    toConsiderEdges[i] = new List<int>();
                toConsiderEdgeSets = new List<int>[size];
                for (int i = 0; i < size; i++)
                    toConsiderEdgeSets[i] = new List<int>();
                toConsiderEdgeColors = new List<int>[size];
                for (int i = 0; i < size; i++)
                    toConsiderEdgeColors[i] = new List<int>();
                toConsiderIntersects = new List<int>[size];
                for (int i = 0; i < size; i++)
                    toConsiderIntersects[i] = new List<int>();
                toConsiderCellColors = new List<int>[size];
                for (int i = 0; i < size; i++)
                    toConsiderCellColors[i] = new List<int>();
            }
            else
            {
                for (int i = 0; i < size; i++)
                {
                    moves[i].Clear();
                    toConsiderEdges[i].Clear();
                    toConsiderEdgeSets[i].Clear();
                    toConsiderEdgeColors[i].Clear();
                    toConsiderIntersects[i].Clear();
                    toConsiderCellColors[i].Clear();
                }
            }

        }

        private bool Perform(IAction? action, List<IAction> backup, int maxDepth, List<int>? initialIntersectsAffected)
        {
            if (maxDepth > edges.Count)
                maxDepth = edges.Count;

            AllocateForPerform(maxDepth + 1);

            if (action != null)
            {
                moves[0].Add(action);
            }

            if (initialIntersectsAffected != null)
                toConsiderIntersects[0].AddRange(initialIntersectsAffected);

            if (edgesSeen == null)
                edgesSeen = new EdgeState[edges.Count];
            else
                Array.Clear(edgesSeen, 0, edgesSeen.Length);
            if (edgePairsSeen == null)
                edgePairsSeen = new TriState[edges.Count, edges.Count];
            else
                ClearEdgePairs();
            if (edgeRestrictsSeen == null)
                edgeRestrictsSeen = new EdgePairRestriction[edges.Count, edges.Count];
            else
                ClearEdgeRestricts();
            if (cellPairsSeen == null)
                cellPairsSeen = new TriState[cells.Count + 1, cells.Count + 1];
            else
                ClearCellPairs();
            if (cellsSeen == null)
                cellsSeen = new bool[cells.Count];
            if (cellColorEdgeColorsSeen == null)
                cellColorEdgeColorsSeen = new bool[cells.Count];
            if (intersectsSeen == null)
                intersectsSeen = new bool[intersections.Count];
            if (interactsSeen == null)
                interactsSeen = new bool[cells.Count];
            if (interactsSeen2 == null)
                interactsSeen2 = new bool[intersections.Count];
            if (interactsSeen3 == null)
                interactsSeen3 = new bool[edges.Count];
            if (colorSetsSeen == null)
                colorSetsSeen = new bool[edges.Count];
            else
            {
                Array.Clear(colorSetsSeen, 0, colorSetsSeen.Length);
            }
            for (int curDepth = 0; curDepth <= maxDepth; curDepth++)
            {
                if (moves[curDepth].Count == 0 &&
                    toConsiderEdges[curDepth].Count == 0 &&
                    toConsiderEdgeSets[curDepth].Count == 0 &&
                    toConsiderEdgeColors[curDepth].Count == 0 &&
                    toConsiderIntersects[curDepth].Count == 0 &&
                    toConsiderCellColors[curDepth].Count == 0
                    )
                    continue;
                if (curDepth > deepestNonEmpty)
                    deepestNonEmpty = curDepth;
                List<IAction> movesCurDepth = moves[curDepth];
                for (var i = 0; i < movesCurDepth.Count; i++)
                {
                    IAction move = movesCurDepth[i];
                    if (!move.Perform())
                        return false;
                    if (superSlowMo)
                    {
                        // Edge restriction actions have no visualization, so no point refreshing.
                        if (!(move is EdgeRestrictionAction))
                        {
                            if (MeshChangeUpdate != null)
                                MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, null, false));
                        }
                    }
                    backup.Add(move);
                    if (!move.Successful)
                        return false;
                }
                if (superSlowMo)
                {
                    if (MeshChangeUpdate != null)
                        MeshChangeUpdate(this, new MeshChangeUpdateEventArgs(this, null, false));
                }
                if (curDepth == maxDepth)
                    break;
                Array.Clear(cellsSeen, 0, cellsSeen.Length);
                Array.Clear(cellColorEdgeColorsSeen, 0, cellColorEdgeColorsSeen.Length);
                Array.Clear(intersectsSeen, 0, intersectsSeen.Length);
                Array.Clear(interactsSeen, 0, interactsSeen.Length);
                Array.Clear(interactsSeen2, 0, interactsSeen2.Length);
                Array.Clear(interactsSeen3, 0, interactsSeen3.Length);

                // Now that moves have been performed, we extract their deepest darkest secrets to work out what needs considering.
                for (var i = 0; i < movesCurDepth.Count; i++)
                {
                    IAction move = movesCurDepth[i];
                    if (move is SetAction)
                    {
                        SetAction setAction = (SetAction) move;
                        AddEdgetToConsider(action, toConsiderEdges, curDepth, setAction.EdgeIndex);
                        var affectedEdges = setAction.GetAffectedEdges();
                        for (var j = 0; j < affectedEdges.Count; j++)
                        {
                            int edge = affectedEdges[j];
                            AddEdgetToConsider(action, toConsiderEdgeSets, curDepth, edge);
                        }
                    }
                    else if (move is ColorJoinAction)
                    {
                        ColorJoinAction cjAction = (ColorJoinAction) move;
                        var affectedEdges = cjAction.GetAffectedEdges();
                        for (var j = 0; j < affectedEdges.Count; j++)
                        {
                            int edge = affectedEdges[j];
                            AddEdgetToConsider(action, toConsiderEdgeColors, curDepth, edge);
                        }
                    }
                    else if (move is CellColorJoinAction)
                    {
                        CellColorJoinAction ccjAction = (CellColorJoinAction) move;
                        var affectedCells = ccjAction.GetAffectedCells();
                        for (var j = 0; j < affectedCells.Count; j++)
                        {
                            int cell = affectedCells[j];
                            AddCellToConsider(action, toConsiderCellColors, curDepth, cell);
                        }
                    }
                    else if (move is EdgeRestrictionAction)
                    {
                        EdgeRestrictionAction erAction = (EdgeRestrictionAction) move;
                        var affectedEdges = erAction.GetAffectedEdges();
                        for (var j = 0; j < affectedEdges.Count; j++)
                        {
                            int edge = affectedEdges[j];
                            AddEdgetToConsider(action, toConsiderEdgeColors, curDepth, edge);
                        }
                    }
                }
                Sort(toConsiderEdges[curDepth]);
                Sort(toConsiderEdgeSets[curDepth]);
                Sort(toConsiderEdgeColors[curDepth]);
                Sort(toConsiderIntersects[curDepth]);
                Sort(toConsiderCellColors[curDepth]);

                if (!ConsiderEdges(moves, toConsiderEdges, curDepth))
                    return false;
                if (!ConsiderIntersects(moves, toConsiderIntersects, curDepth))
                    return false;
                if (!ConsiderEdgeSets(moves, toConsiderEdgeSets, curDepth))
                    return false;
                if (!ConsiderEdgeColors(moves, toConsiderEdgeColors, curDepth))
                    return false;
                if (!ConsiderCellColors(moves, toConsiderCellColors, curDepth))
                    return false;

            }
            return true;
        }
        private void Sort(List<int> list)
        {
            list.Sort();
        }
        private void ClearCellPairs()
        {
            for (var i = 0; i < cellPairsToClean.Count; i++)
            {
                KeyValuePair<int, int> kvp = cellPairsToClean[i];
                cellPairsSeen[kvp.Key, kvp.Value] = TriState.Unknown;
                cellPairsSeen[kvp.Value, kvp.Key] = TriState.Unknown;
            }
            cellPairsToClean.Clear();
        }

        private void ClearEdgeRestricts()
        {
            for (var i = 0; i < edgeRestrictsToClean.Count; i++)
            {
                KeyValuePair<int, int> kvp = edgeRestrictsToClean[i];
                edgeRestrictsSeen[kvp.Key, kvp.Value] = EdgePairRestriction.None;
                edgeRestrictsSeen[kvp.Value, kvp.Key] = EdgePairRestriction.None;
            }
            edgeRestrictsToClean.Clear();
        }

        private void ClearEdgePairs()
        {
            for (var i = 0; i < edgePairsToClean.Count; i++)
            {
                KeyValuePair<int, int> kvp = edgePairsToClean[i];
                edgePairsSeen[kvp.Key, kvp.Value] = TriState.Unknown;
                edgePairsSeen[kvp.Value, kvp.Key] = TriState.Unknown;
            }
            edgePairsToClean.Clear();
        }

        private bool ConsiderCellColors(List<IAction>[] moves, List<int>[] toConsiderCellColors, int curDepth)
        {
            if (useCellColoring)
            {
                int lastCell = -1;
                for (var i = 0; i < toConsiderCellColors[curDepth].Count; i++)
                {
                    int cellColorIndex = toConsiderCellColors[curDepth][i];
                    if (cellColorIndex == lastCell)
                        continue;
                    lastCell = cellColorIndex;
                    Cell cell = cells[cellColorIndex];
                    if (!GatherCellColoringMoves(cell, moves, curDepth, edgesSeen, cellColorIndex))
                        return false;
                    // TODO: can cell color interact with intersection types?  Edge coloring interactions will probably cover it, so long as they are turned on.
                    if (!cellColorEdgeColorsSeen[cellColorIndex])
                    {
                        cellColorEdgeColorsSeen[cellColorIndex] = true;
                        if (!GatherCellColoringEdgeColoringMovesForCellColorChange(cell, moves, curDepth,
                            cellColorIndex))
                            return false;
                    }
                }
            }
            return true;
        }

        private bool ConsiderEdgeColors(List<IAction>[] moves, List<int>[] toConsiderEdgeColors, int curDepth)
        {
            if (UseColoring || UseEdgeRestricts)
            {
                int lastEdge = -1;
                for (var i = 0; i < toConsiderEdgeColors[curDepth].Count; i++)
                {
                    int edgeAffectedIndex = toConsiderEdgeColors[curDepth][i];
                    if (edgeAffectedIndex == lastEdge)
                        continue;
                    lastEdge = edgeAffectedIndex;
                    Edge edge = edges[edgeAffectedIndex];
                    if (UseColoring && edge.Color != 0)
                    {
                        // TODO: do not activate this unless the edge color has changed.
                        GatherCellColoringEdgeColoringMovesForEdgeColorChange(edge, moves, curDepth, edgeAffectedIndex);
                    }
                    foreach (int intersIndex in edge.Intersections)
                    {
                        if (interactsSeen2[intersIndex])
                            continue;
                        else
                            interactsSeen2[intersIndex] = true;
                        Intersection inters = intersections[intersIndex];
                        if (!GatherInteractForcedMoves(inters, moves, curDepth, edgesSeen, edgePairsSeen,
                            edgeRestrictsSeen))
                            return false;
                    }
                    if (edge.Color != 0)
                    {
                        if (!GatherFollowColoringColorSetChanged(moves, curDepth, edge, edgesSeen, edgeAffectedIndex))
                            return false;
                    }
                    if (!GatherFollowEdgeRestrictions(moves, curDepth, edge, edgesSeen, edgeAffectedIndex))
                        return false;
                    /*
                    if (UseCellPairs || UseCellPairsTopLevel)
                    {
                        foreach (int intersIndex in edge.Intersections)
                        {
                            Intersection inters = intersections[intersIndex];
                            foreach (int divider in inters.Edges)
                            {
                                if (!interactsSeen3[divider])
                                {
                                    interactsSeen3[divider] = true;
                                    Edge dividingEdge = edges[divider];
                                    if (!GatherCellPairForcedMoves(dividingEdge, moves, curDepth, edgesSeen, edgePairsSeen, edgeRestrictsSeen))
                                        return false;
                                }
                            }
                        }
                        foreach (int cellIndex in edge.Cells)
                        {
                            Cell cell = cells[cellIndex];
                            foreach (int divider in cell.Edges)
                            {
                                if (!interactsSeen3[divider])
                                {
                                    interactsSeen3[divider] = true;
                                    Edge dividingEdge = edges[divider];
                                    if (!GatherCellPairForcedMoves(dividingEdge, moves, curDepth, edgesSeen, edgePairsSeen, edgeRestrictsSeen))
                                        return false;
                                }
                            }
                        }
                    }
                     */
                }
            }
            return true;
        }

        private bool ConsiderEdgeSets(List<IAction>[] moves, List<int>[] toConsiderEdgeSets, int curDepth)
        {
            if (considerMultipleLoops)
            {
                int lastEdge = -1;
                for (var i = 0; i < toConsiderEdgeSets[curDepth].Count; i++)
                {
                    int edgeAffectedIndex = toConsiderEdgeSets[curDepth][i];
                    if (edgeAffectedIndex == lastEdge)
                        continue;
                    lastEdge = edgeAffectedIndex;
                    Edge edge = edges[edgeAffectedIndex];
                    if (!GatherExcludeClosingLoopEarly(moves, curDepth, edgeAffectedIndex, edge, edgesSeen))
                        return false;
                }
            }
            return true;
        }

        private bool ConsiderIntersects(List<IAction>[] moves, List<int>[] toConsiderIntersects, int curDepth)
        {
            int lastIntersect = -1;
            for (var i = 0; i < toConsiderIntersects[curDepth].Count; i++)
            {
                int intersectIndex = toConsiderIntersects[curDepth][i];
                if (intersectIndex == lastIntersect)
                    continue;
                lastIntersect = intersectIndex;
                Intersection inters = intersections[intersectIndex];
                if (!intersectsSeen[intersectIndex])
                {
                    intersectsSeen[intersectIndex] = true;
                    if (!GatherIntersectionForcedMoves(inters, moves, curDepth, edgesSeen, intersectIndex))
                        return false;
                }

                if (UseColoring || UseEdgeRestricts)
                {
                    if (!interactsSeen2[intersectIndex])
                    {
                        interactsSeen2[intersectIndex] = true;
                        if (!GatherInteractForcedMoves(inters, moves, curDepth, edgesSeen, edgePairsSeen,
                            edgeRestrictsSeen))
                            return false;
                    }
                }
                /*
                if (UseCellPairs || UseCellPairsTopLevel)
                {
                    foreach (int divider in inters.Edges)
                    {
                        if (!interactsSeen3[divider])
                        {
                            interactsSeen3[divider] = true;
                            Edge dividingEdge = edges[divider];
                            if (!GatherCellPairForcedMoves(dividingEdge, moves, curDepth, edgesSeen, edgePairsSeen, edgeRestrictsSeen))
                                return false;
                        }
                    }
                }
                 */
            }
            return true;
        }

        private bool ConsiderEdges(List<IAction>[] moves, List<int>[] toConsiderEdges, int curDepth)
        {
            int lastEdge = -1;
            for (var i = 0; i < toConsiderEdges[curDepth].Count; i++)
            {
                int edgeAffectedIndex = toConsiderEdges[curDepth][i];
                if (edgeAffectedIndex == lastEdge)
                    continue;
                lastEdge = edgeAffectedIndex;

                Edge edge = edges[edgeAffectedIndex];

                for (var j = 0; j < edge.Intersections.Length; j++)
                {
                    int intersIndex = edge.Intersections[j];
                    if (intersectsSeen[intersIndex])
                        continue;
                    else
                        intersectsSeen[intersIndex] = true;
                    Intersection inters = intersections[intersIndex];
                    if (!GatherIntersectionForcedMoves(inters, moves, curDepth, edgesSeen, intersIndex))
                        return false;
                }
                if (UseColoring || UseEdgeRestricts)
                {
                    for (var j = 0; j < edge.Intersections.Length; j++)
                    {
                        int intersIndex = edge.Intersections[j];
                        if (interactsSeen2[intersIndex])
                            continue;
                        else
                            interactsSeen2[intersIndex] = true;
                        Intersection inters = intersections[intersIndex];
                        if (!GatherInteractForcedMoves(inters, moves, curDepth, edgesSeen, edgePairsSeen,
                            edgeRestrictsSeen))
                            return false;
                    }
                }
                /*
                if (UseCellPairs || UseCellPairsTopLevel)
                {
                    foreach (int intersIndex in edge.Intersections)
                    {
                        Intersection inters = intersections[intersIndex];
                        foreach (int divider in inters.Edges)
                        {
                            if (!interactsSeen3[divider])
                            {
                                interactsSeen3[divider] = true;
                                Edge dividingEdge = edges[divider];
                                if (!GatherCellPairForcedMoves(dividingEdge, moves, curDepth, edgesSeen, edgePairsSeen, edgeRestrictsSeen))
                                    return false;
                            }
                        }
                    }
                    foreach (int cellIndex in edge.Cells)
                    {
                        Cell cell = cells[cellIndex];
                        foreach (int divider in cell.Edges)
                        {
                            if (!interactsSeen3[divider])
                            {
                                interactsSeen3[divider] = true;
                                Edge dividingEdge = edges[divider];
                                if (!GatherCellPairForcedMoves(dividingEdge, moves, curDepth, edgesSeen, edgePairsSeen, edgeRestrictsSeen))
                                    return false;
                            }
                        }
                    }
                }
                 */
                if (!GatherCellColoringMoves(edge, moves, curDepth, edgesSeen, edgeAffectedIndex))
                    return false;
                if (UseColoring)
                {
                    if (edge.Color != 0)
                    {
                        if (!GatherFollowColoring(moves, curDepth, edge, edgesSeen, edgeAffectedIndex))
                            return false;
                    }
                }
                if (!GatherFollowEdgeRestrictions(moves, curDepth, edge, edgesSeen, edgeAffectedIndex))
                    return false;
            }
            return true;
        }

        private void AddEdgetToConsider(IAction? sourceAction, List<int>[] toConsider, int curDepth, int other)
        {
            int depth = curDepth;
            if (sourceAction != null)
            {
                int dist = GetEdgeDistance(other, sourceAction);
                if (dist > curDepth)
                    depth = dist;
            }
            if (depth < toConsider.Length)
                toConsider[depth].Add(other);
        }

        private void AddCellToConsider(IAction? sourceAction, List<int>[] toConsider, int curDepth, int cell)
        {
            int depth = curDepth;
            if (sourceAction != null)
            {
                int minDist = GetCellDistance(cell, sourceAction);
                if (minDist > curDepth)
                    depth = minDist;
            }
            if (depth < toConsider.Length)
                toConsider[depth].Add(cell);
        }

        private bool GatherExcludeClosingLoopEarly(List<IAction>[] moves, int curDepth, int edgeAffectedIndex, Edge edge, EdgeState[] edgesSeen)
        {
            for (var i2 = 0; i2 < edge.Intersections.Length; i2++)
            {
                int intersIndex = edge.Intersections[i2];
                Intersection inters = intersections[intersIndex];
                if (inters.FilledCount != 1)
                    continue;
                for (var j = 0; j < inters.Edges.Count; j++)
                {
                    int otherEdgeIndex = inters.Edges[j];
// If this edge is still empty, there is no reason why we shouldn't check it for closing the loop early.
                    //if (otherEdgeIndex == edgeAffectedIndex)
                    //    continue;
                    Edge otherEdge = edges[otherEdgeIndex];
                    if (otherEdge.State == EdgeState.Empty)
                    {
                        int edgeSet1 = GetEdgeSet(otherEdge.Intersections[0], otherEdgeIndex);
                        int edgeSet2 = GetEdgeSet(otherEdge.Intersections[1], otherEdgeIndex);
                        if (edgeSet1 != 0 && edgeSet1 == edgeSet2)
                        {
                            // Would close a loop if we did it.
                            bool okay = false;
                            int nonEmpty = 0;
                            foreach (List<int> edgeSet in edgeSets)
                                if (edgeSet.Count > 0)
                                    nonEmpty++;
                            if (nonEmpty == 1)
                            {
                                Intersection inters1 = intersections[otherEdge.Intersections[0]];
                                Intersection inters2 = intersections[otherEdge.Intersections[1]];
                                if (inters1.FilledCount == 1 && inters2.FilledCount == 1)
                                {
                                    if (satisifiedIntersCount == intersections.Count - 2)
                                    {
                                        // TODO: make this less expensive.
                                        bool fine = true;
                                        for (int i = 0; i < intersections.Count; i++)
                                        {
                                            if (intersections[i].Type != IntersType.Unknown)
                                            {
                                                if (DetermineTypeAssuming(i, otherEdgeIndex) != intersections[i].Type)
                                                {
                                                    fine = false;
                                                    break;
                                                }
                                            }
                                        }
                                        if (fine)
                                            okay = true;
                                    }
                                }
                            }
                            if (!okay)
                            {
                                if (edgesSeen[otherEdgeIndex] == EdgeState.Filled)
                                    return false;
                                if (edgesSeen[otherEdgeIndex] == EdgeState.Excluded)
                                    continue;
                                edgesSeen[otherEdgeIndex] = EdgeState.Excluded;
                                moves[curDepth + 1].Add(new SetAction(this, otherEdgeIndex, EdgeState.Excluded));
                            }
                        }
                    }
                }
            }
            return true;
        }

        private bool GatherFollowColoringColorSetChanged(List<IAction>[] moves, int curDepth, Edge edge, EdgeState[] edgesSeen, int edgeIndex)
        {
            List<int> colorSet;
            colorSet = colorSets[Math.Abs(edge.Color) - 1];
            // Stop caching that we have seen this if we were.
            colorSetsSeen[Math.Abs(edge.Color) - 1] = false;
            EdgeState posState = EdgeState.Empty;
            for (var index = 0; index < colorSet.Count; index++)
            {
                int i = colorSet[index];
                if (i == edgeIndex)
                    continue;
                Edge toCheck = edges[i];
                if (toCheck.State != EdgeState.Empty)
                {
                    if (toCheck.State == EdgeState.Filled)
                        posState = toCheck.Color > 0 ? EdgeState.Filled : EdgeState.Excluded;
                    else
                        posState = toCheck.Color > 0 ? EdgeState.Excluded : EdgeState.Filled;
                    break;
                }
            }
            if (posState != EdgeState.Empty)
            {
                EdgeState expectedState = edge.Color > 0 ? posState : (posState == EdgeState.Excluded ? EdgeState.Filled : EdgeState.Excluded);
                if (edge.State == EdgeState.Empty)
                {
                    if (!AddSetAction(edgeIndex, expectedState, moves, curDepth + 1, edgesSeen))
                        return false;
                }
                else if (edge.State != expectedState)
                {
                    return false;
                }
            }
            return true;
        }

        private bool GatherFollowColoring(List<IAction>[] moves, int curDepth, Edge edge, EdgeState[] edgesSeen, int edgeIndex)
        {
            List<int> colorSet;
            colorSet = colorSets[Math.Abs(edge.Color) - 1];
            if (colorSetsSeen[Math.Abs(edge.Color) - 1])
                return true;
            colorSetsSeen[Math.Abs(edge.Color) - 1] = true;
            for (var index = 0; index < colorSet.Count; index++)
            {
                int i = colorSet[index];
                int dist = edgeDistances[i, edgeIndex];
                if (curDepth + dist >= moves.Length)
                    continue;
                Edge toCheck = edges[i];
                if (toCheck.Color == edge.Color && toCheck.State != edge.State)
                {
                    if (toCheck.State != EdgeState.Empty)
                        return false;
                    EdgeState newState = edge.State;
                    if (!AddSetAction(i, newState, moves, curDepth + dist, edgesSeen))
                        return false;
                }
                else if (toCheck.Color == -edge.Color)
                {
                    if (toCheck.State == edge.State)
                        return false;
                    if (toCheck.State != EdgeState.Empty)
                        continue;
                    EdgeState newState = edge.State == EdgeState.Filled ? EdgeState.Excluded : EdgeState.Filled;
                    if (!AddSetAction(i, newState, moves, curDepth + dist, edgesSeen))
                        return false;
                }
            }
            return true;
        }

        private bool GatherFollowEdgeRestrictions(List<IAction>[] moves, int curDepth, Edge edge, EdgeState[] edgesSeen, int edgeIndex)
        {
            if (!UseEdgeRestricts || edgePairRestrictions == null)
                return true;
            if (edge.State == EdgeState.Empty)
                return true;
            for (var i = 0; i < edge.Intersections.Length; i++)
            {
                int inters = edge.Intersections[i];
                Intersection inter = intersections[inters];
                for (var j = 0; j < inter.Edges.Count; j++)
                {
                    int otherEdgeIndex = inter.Edges[j];
                    if (otherEdgeIndex == edgeIndex)
                        continue;
                    Edge toCheck = edges[otherEdgeIndex];
                    switch (edgePairRestrictions[edgeIndex, otherEdgeIndex])
                    {
                        case EdgePairRestriction.NotBoth:
                            if (edge.State == EdgeState.Filled)
                            {
                                if (toCheck.State == EdgeState.Filled)
                                    return false;
                                if (toCheck.State != EdgeState.Empty)
                                    continue;
                                if (!AddSetAction(otherEdgeIndex, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                            break;
                        case EdgePairRestriction.NotNeither:
                            if (edge.State == EdgeState.Excluded)
                            {
                                if (toCheck.State == EdgeState.Excluded)
                                    return false;
                                if (toCheck.State != EdgeState.Empty)
                                    continue;
                                if (!AddSetAction(otherEdgeIndex, EdgeState.Filled, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                            break;
                    }
                }
            }
            for (var i = 0; i < edge.Cells.Count; i++)
            {
                int cellIndex = edge.Cells[i];
                Cell cell = cells[cellIndex];
                for (var j = 0; j < cell.Edges.Count; j++)
                {
                    int otherEdgeIndex = cell.Edges[j];
                    if (otherEdgeIndex == edgeIndex)
                        continue;
                    Edge toCheck = edges[otherEdgeIndex];
                    switch (edgePairRestrictions[edgeIndex, otherEdgeIndex])
                    {
                        case EdgePairRestriction.NotBoth:
                            if (edge.State == EdgeState.Filled)
                            {
                                if (toCheck.State == EdgeState.Filled)
                                    return false;
                                if (toCheck.State != EdgeState.Empty)
                                    continue;
                                if (!AddSetAction(otherEdgeIndex, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                            break;
                        case EdgePairRestriction.NotNeither:
                            if (edge.State == EdgeState.Excluded)
                            {
                                if (toCheck.State == EdgeState.Excluded)
                                    return false;
                                if (toCheck.State != EdgeState.Empty)
                                    continue;
                                if (!AddSetAction(otherEdgeIndex, EdgeState.Filled, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                            break;
                    }
                }
            }
            return true;
        }

        private bool AddSetAction(int edgeIndex, EdgeState newState, List<IAction>[] moves, int targetDepth, EdgeState[] edgesSeen)
        {
            if (edgesSeen[edgeIndex] != EdgeState.Empty && edgesSeen[edgeIndex] != newState)
                return false;
            if (edgesSeen[edgeIndex] == newState)
                return true;
            edgesSeen[edgeIndex] = newState;
            moves[targetDepth].Add(new SetAction(this, edgeIndex, newState));
            return true;
        }

        private bool AddColorJoinAction(int edge1, int edge2, bool same, List<IAction>[] moves, int targetDepth, TriState[,] colorJoinsSeen)
        {
            TriState newTriState = same ? TriState.Same : TriState.Opposite;
            if (colorJoinsSeen[edge1, edge2] != TriState.Unknown && colorJoinsSeen[edge1, edge2] != newTriState)
                return false;
            if (colorJoinsSeen[edge1, edge2] == newTriState)
                return true;
            colorJoinsSeen[edge1, edge2] = newTriState;
            colorJoinsSeen[edge2, edge1] = newTriState;
            edgePairsToClean.Add(new KeyValuePair<int, int>(edge1, edge2));
            moves[targetDepth].Add(new ColorJoinAction(this, edge1, edge2, same));
            return true;
        }

        private bool AddEdgeRestrictAction(int edge1, int edge2, EdgePairRestriction edgePairRestriction, List<IAction>[] moves, int targetDepth, EdgePairRestriction[,] edgeRestrictsSeen, TriState[,] edgePairsSeen)
        {
            if (edgeRestrictsSeen[edge1, edge2] != EdgePairRestriction.None && edgeRestrictsSeen[edge1, edge2] != edgePairRestriction)
            {
                // Edge restrictions are special, seing both options isn't a contradiction.
                // It means we've discovered a color info.
                // return AddColorJoinAction(edge1, edge2, false, moves, targetDepth, edgePairsSeen);
                // But lets let it be discovered by the color code instead.
                return true;
            }
            if (edgeRestrictsSeen[edge1, edge2] == edgePairRestriction)
                return true;
            // TODO: remove this once the logic no longer produces results we've already got.
            if (edgePairRestrictions != null && edgePairRestrictions[edge1, edge2] != EdgePairRestriction.None)
                return true;
            edgeRestrictsSeen[edge1, edge2] = edgePairRestriction;
            edgeRestrictsSeen[edge2, edge1] = edgePairRestriction;
            edgeRestrictsToClean.Add(new KeyValuePair<int, int>(edge1, edge2));
            moves[targetDepth].Add(new EdgeRestrictionAction(this, edge1, edge2, edgePairRestriction));
            return true;
        }

        private bool AddCellColorJoinAction(int cell1, int cell2, bool same, List<IAction>[] moves, int targetDepth, TriState[,] cellColorJoinsSeen)
        {
            if (cell1 == -1)
            {
                int temp = cell2;
                cell2 = cell1;
                cell1 = temp;
            }
            TriState newTriState = same ? TriState.Same : TriState.Opposite;
            if (cellColorJoinsSeen[cell1 + 1, cell2 + 1] != TriState.Unknown && cellColorJoinsSeen[cell1 + 1, cell2 + 1] != newTriState)
                return false;
            if (cellColorJoinsSeen[cell1 + 1, cell2 + 1] == newTriState)
                return true;
            cellColorJoinsSeen[cell1 + 1, cell2 + 1] = newTriState;
            cellColorJoinsSeen[cell2 + 1, cell1 + 1] = newTriState;
            cellPairsToClean.Add(new KeyValuePair<int, int>(cell1 + 1, cell2 + 1));
            moves[targetDepth].Add(new CellColorJoinAction(this, cell1, cell2, same));
            return true;
        }

        private bool GatherIntersectionForcedMoves(Intersection inters, List<IAction>[] moves, int curDepth, EdgeState[] edgesSeen, int intersIndex)
        {
            bool definiteLines = false;
            if (inters.Type != IntersType.Unknown)
                definiteLines = true;
            EdgeState toPerform = EdgeState.Empty;
            if (inters.FilledCount == 2 && inters.ExcludedCount < inters.Edges.Count - 2 || inters.FilledCount == 0 && inters.ExcludedCount > inters.Edges.Count - 2 && inters.ExcludedCount < inters.Edges.Count)
            {
                toPerform = EdgeState.Excluded;
            }
            if (inters.Edges.Count - inters.ExcludedCount == 2 && inters.FilledCount < 2 && (inters.FilledCount > 0 || definiteLines))
            {
                toPerform = EdgeState.Filled;
            }
            if (toPerform != EdgeState.Empty)
            {
                for (var i = 0; i < inters.Edges.Count; i++)
                {
                    int otherEdgeIndex = inters.Edges[i];
                    Edge otherEdge = edges[otherEdgeIndex];
                    if (otherEdge.State == EdgeState.Empty)
                    {
                        if (!AddSetAction(otherEdgeIndex, toPerform, moves, curDepth + 1, edgesSeen))
                            return false;
                    }
                }
            }

            // Special case of IntersType.Turn, blank intersection next to it is interesting when one of its edges becomes excluded.
            if (inters.ExcludedCount > 0 && inters.ExcludedCount < inters.Edges.Count)
            {
                for (int i = 0; i < inters.Edges.Count; i++)
                {
                    int edge = inters.Edges[i];
                    if (edges[edge].State == EdgeState.Excluded)
                    {
                        int oppositeEdge = inters.Opposite[i];
                        if (oppositeEdge != -1 && edges[oppositeEdge].State != EdgeState.Excluded)
                        {
                            int otherInters = GetOtherInters(oppositeEdge, intersIndex);
                            Intersection otherIntersection = intersections[otherInters];
                            if (otherIntersection.Type == IntersType.Turn)
                            {
                                if (edges[oppositeEdge].State == EdgeState.Filled)
                                    return false;
                                if (!AddSetAction(oppositeEdge, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                        }
                    }
                }
            }
            // Special case of IntersType.Turn, blank intersection next to it is interesting when one of its perpendicular edges becomes filled.
            if (inters.FilledCount == 1)
            {
                for (int i = 0; i < inters.Edges.Count; i++)
                {
                    int edge = inters.Edges[i];
                    if (edges[edge].State == EdgeState.Filled)
                    {
                        for (int j = 0; j < inters.Edges.Count; j++)
                        {
                            if (j == i)
                                continue;
                            if (inters.Opposite[j] == edge)
                                continue;
                            int perpEdge = inters.Edges[j];
                            if (edges[perpEdge].State != EdgeState.Excluded)
                            {
                                int otherInters = GetOtherInters(perpEdge, intersIndex);
                                Intersection otherIntersection = intersections[otherInters];
                                if (otherIntersection.Type == IntersType.Turn)
                                {
                                    // This won't happen...
                                    if (edges[perpEdge].State == EdgeState.Filled)
                                        return false;
                                    // But this might :)
                                    if (!AddSetAction(perpEdge, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                        return false;
                                }
                            }
                        }
                    }
                }
            }
                // Special case of InterType.Twist, blank intersection next to it is interesting when the straight through line is extended through it.
            else if (inters.FilledCount == 2)
            {
                for (int i = 0; i < inters.Edges.Count; i++)
                {
                    int edge = inters.Edges[i];
                    if (edges[edge].State == EdgeState.Filled)
                    {
                        int opposite = inters.Opposite[i];
                        if (opposite != -1 && edges[opposite].State == EdgeState.Filled)
                        {
                            int otherInters = GetOtherInters(opposite, intersIndex);
                            Intersection otherIntersection = intersections[otherInters];
                            if (otherIntersection.Type == IntersType.Twist)
                            {
                                // Onwards ... we have so far to go...
                                for (int j = 0; j < otherIntersection.Edges.Count; j++)
                                {
                                    if (otherIntersection.Edges[j] == opposite)
                                    {
                                        int moreOpposite = otherIntersection.Opposite[j];
                                        if (moreOpposite != -1)
                                        {
                                            int thirdInters = GetOtherInters(moreOpposite, otherInters);
                                            Intersection thirdIntersection = intersections[thirdInters];
                                            for (int k = 0; k < thirdIntersection.Edges.Count; k++)
                                            {
                                                if (thirdIntersection.Edges[k] == moreOpposite)
                                                {
                                                    int thirdOp = thirdIntersection.Opposite[k];
                                                    if (thirdOp != -1 && edges[thirdOp].State != EdgeState.Excluded)
                                                    {
                                                        if (edges[thirdOp].State == EdgeState.Filled)
                                                            return false;
                                                        if (!AddSetAction(thirdOp, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                                            return false;
                                                    }
                                                }
                                            }
                                        }
                                        else
                                            // This will surely have been handled elsewhere, but ... anyway.
                                            return false;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (inters.Type == IntersType.Turn)
            {
                if (inters.FilledCount == 1)
                {
                    for (int i=0; i < inters.Edges.Count; i++)
                    {
                        int edge = inters.Edges[i];
                        if (edges[edge].State == EdgeState.Filled)
                        {
                            int other = inters.Opposite[i];
                            if (other != -1 && edges[other].State != EdgeState.Excluded)
                            {
                                if (edges[other].State == EdgeState.Filled)
                                    return false;
                                if (!AddSetAction(other, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                        }
                    }
                }
                if ((inters.ExcludedCount == 1 && inters.Edges.Count == 4) || inters.Edges.Count == 3)
                {
                    for (int i = 0; i < inters.Edges.Count; i++)
                    {
                        int opposite = inters.Opposite[i];
                        int edge = inters.Edges[i];
                        if (opposite == -1 || edges[opposite].State == EdgeState.Excluded)
                        {
                            if (edges[edge].State != EdgeState.Filled)
                            {
                                if (edges[edge].State == EdgeState.Excluded)
                                    return false;
                                if (!AddSetAction(edge, EdgeState.Filled, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                        }
                    }
                }
                
                // TODO: this code path is static - as the excluded edge case is actually handled seperately, above.  We can optimize this to be run at top level only?
                for (int i = 0; i < inters.Edges.Count; i++)
                {
                    int edge = inters.Edges[i];
                    int otherInters = GetOtherInters(edge, intersIndex);
                    Intersection otherIntersection = intersections[otherInters];
                    for (int j = 0; j < otherIntersection.Edges.Count; j++)
                    {
                        if (otherIntersection.Edges[j] == edge)
                        {
                            int other = otherIntersection.Opposite[j];
                            if (other == -1 || edges[other].State == EdgeState.Excluded)
                            {
                                if (edges[edge].State != EdgeState.Excluded)
                                {
                                    if (edges[edge].State == EdgeState.Filled)
                                        return false;
                                    if (!AddSetAction(edge, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                        return false;
                                }
                            }
                            // TODO: in the else case of the above we should add an edge restriction, but it is not of a supported type.  It is of the type NotFilledExcluded, and when combined with 
                            // edge restriction processing at that intersection would be very useful.
                        }
                    }
                }
                if (UseColoring)
                {
                    // TODO: add anticoloring between opposite edges - not really needed as it will be derived by trivial trials if DeriveColoring is on.
                }
                if (inters.FilledCount > 0)
                {
                    for (int i = 0; i < inters.Edges.Count; i++)
                    {
                        int edge = inters.Edges[i];
                        if (edges[edge].State == EdgeState.Filled)
                        {
                            int otherInter = GetOtherInters(edge, intersIndex);
                            Intersection otherIntersection = intersections[otherInter];
                            for (int j = 0; j < otherIntersection.Edges.Count; j++)
                            {
                                if (otherIntersection.Edges[j] == edge)
                                {
                                    int toBeSet = otherIntersection.Opposite[j];
                                    if (toBeSet == -1)
                                        return false;
                                    if (edges[toBeSet].State == EdgeState.Excluded)
                                        return false;
                                    if (edges[toBeSet].State != EdgeState.Filled)
                                    {
                                        if (!AddSetAction(toBeSet, EdgeState.Filled, moves, curDepth + 1, edgesSeen))
                                            return false;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if (inters.Type == IntersType.Twist)
            {
                if (inters.FilledCount == 1)
                {
                    for (int i = 0; i < inters.Edges.Count; i++)
                    {
                        int edge = inters.Edges[i];
                        if (edges[edge].State == EdgeState.Filled)
                        {
                            int other = inters.Opposite[i];
                            if (other != -1 && edges[other].State != EdgeState.Filled)
                            {
                                if (edges[other].State == EdgeState.Excluded)
                                    return false;
                                if (!AddSetAction(other, EdgeState.Filled, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                        }
                    }
                }
                if ((inters.ExcludedCount == 1 && inters.Edges.Count == 4) || inters.Edges.Count == 3)
                {
                    for (int i = 0; i < inters.Edges.Count; i++)
                    {
                        int opposite = inters.Opposite[i];
                        int edge = inters.Edges[i];
                        if (opposite == -1 || edges[opposite].State == EdgeState.Excluded)
                        {
                            if (edges[edge].State != EdgeState.Excluded)
                            {
                                if (edges[edge].State == EdgeState.Filled)
                                    return false;
                                if (!AddSetAction(edge, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                        }
                    }
                }
                if (UseColoring)
                {
                    // TODO: ensure the coloring restrictions are in place as implied for twists.
                    // Not really needed if derive colouring is on, trival trials will discover the associated coloring.
                }
                // If we find a straight, the opposite must be bent.
                if (inters.FilledCount > 0)
                {
                    for (int i = 0; i < inters.Edges.Count; i++)
                    {
                        int edge = inters.Edges[i];
                        if (edges[edge].State == EdgeState.Filled)
                        {
                            int otherInter = GetOtherInters(edge, intersIndex);
                            Intersection otherIntersection = intersections[otherInter];
                            for (int j = 0; j < otherIntersection.Edges.Count; j++)
                            {
                                if (otherIntersection.Edges[j] == edge)
                                {
                                    int checkEdge = otherIntersection.Opposite[j];
                                    if (checkEdge != -1 && edges[checkEdge].State == EdgeState.Filled)
                                    {
                                        // We have a straight, the other one then must be a turn.
                                        int opposite = inters.Opposite[i];
                                        // This should already have been detected, but moving along...
                                        if (opposite == -1)
                                            return false;
                                        int finalInter = GetOtherInters(opposite, intersIndex);
                                        Intersection finalIntersection = intersections[finalInter];
                                        for (int k = 0; k < finalIntersection.Edges.Count; k++)
                                        {
                                            if (finalIntersection.Edges[k] == opposite)
                                            {
                                                int toBeCleared = finalIntersection.Opposite[k];
                                                if (toBeCleared == -1)
                                                    continue;
                                                if (edges[toBeCleared].State == EdgeState.Filled)
                                                    return false;
                                                if (edges[toBeCleared].State != EdgeState.Excluded)
                                                {
                                                    if (!AddSetAction(toBeCleared, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                                        return false;
                                                }

                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return true;
        }

        private bool GatherCellColoringMoves(Cell cell, List<IAction>[] moves, int curDepth, EdgeState[] edgesSeen, int cellIndex)
        {
            if (useCellColoring)
            {
                for (var i = 0; i < cell.Edges.Count; i++)
                {
                    int edge = cell.Edges[i];
                    Edge e = edges[edge];
                    if (!GatherCellColoringMoves(e, moves, curDepth, edgesSeen, edge))
                        return false;
                }
            }
            return true;
        }

        private bool GatherCellColoringMoves(Edge e, List<IAction>[] moves, int curDepth, EdgeState[] edgesSeen, int edgeIndex)
        {
            if (useCellColoring)
            {
                Cell? cell = null;
                int cellIndex = -2;
                Cell? otherCell = null;
                int otherC = -2;
                for (var i = 0; i < e.Cells.Count; i++)
                {
                    int c = e.Cells[i];
                    if (cellIndex == -2)
                    {
                        cellIndex = c;
                        cell = cells[cellIndex];
                    }
                    else
                    {
                        otherC = c;
                        otherCell = cells[otherC];
                    }
                }
                if (cell == null)
                    return true;
                if (otherC != -2 && otherCell != null)
                {
                    if (e.State == EdgeState.Empty)
                    {
                        if (cell.Color != 0)
                        {
                            if (otherCell.Color == cell.Color)
                            {
                                if (!AddSetAction(edgeIndex, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                            else if (otherCell.Color == -cell.Color)
                            {
                                if (!AddSetAction(edgeIndex, EdgeState.Filled, moves, curDepth + 1, edgesSeen))
                                    return false;
                            }
                        }
                    }
                    else
                    {
                        if (cell.Color != 0)
                        {
                            if (otherCell.Color == cell.Color)
                            {
                                if (e.State == EdgeState.Filled)
                                    return false;
                            }
                            else if (otherCell.Color == -cell.Color)
                            {
                                if (e.State == EdgeState.Excluded)
                                    return false;
                            }
                            else
                            {
                                if (!AddCellColorJoinAction(cellIndex, otherC, e.State == EdgeState.Excluded, moves, curDepth + 1, cellPairsSeen))
                                    return false;
                            }
                        }
                        else
                        {
                            if (!AddCellColorJoinAction(cellIndex, otherC, e.State == EdgeState.Excluded, moves, curDepth + 1, cellPairsSeen))
                                return false;
                        }
                    }
                }
                else
                {
                    if (e.State == EdgeState.Empty)
                    {
                        if (cell.Color == 1)
                        {
                            if (!AddSetAction(edgeIndex, EdgeState.Excluded, moves, curDepth + 1, edgesSeen))
                                return false;
                        }
                        else if (cell.Color == -1)
                        {
                            if (!AddSetAction(edgeIndex, EdgeState.Filled, moves, curDepth + 1, edgesSeen))
                                return false;
                        }
                    }
                    else
                    {
                        if (cell.Color != 0)
                        {
                            if (cell.Color == 1)
                            {
                                if (e.State == EdgeState.Filled)
                                    return false;
                            }
                            else if (cell.Color == -1)
                            {
                                if (e.State == EdgeState.Excluded)
                                    return false;
                            }
                            else
                            {
                                if (!AddCellColorJoinAction(cellIndex, -1, e.State == EdgeState.Excluded, moves, curDepth + 1, cellPairsSeen))
                                    return false;
                            }
                        }
                        else
                        {
                            if (!AddCellColorJoinAction(cellIndex, -1, e.State == EdgeState.Excluded, moves, curDepth + 1, cellPairsSeen))
                                return false;
                        }
                    }
                }
            }
            return true;
        }


        private int GetAdjacentCell(int cell, int byEdge)
        {
            List<int> candidates = edges[byEdge].Cells;
            for (int i = 0; i < candidates.Count; i++)
            {
                if (candidates[i] != cell)
                    return candidates[i];
            }
            return -1;
        }

        private bool GatherCellColoringEdgeColoringMovesForCellColorChange(Cell cell, List<IAction>[] moves, int curDepth, int cellIndex)
        {
            if (UseCellColoring && UseColoring)
            {
                for (int i2 = 0; i2 < cell.Edges.Count; i2++)
                {
                    int edge = cell.Edges[i2];
                    Edge e = edges[edge];
                    int otherC = GetAdjacentCell(cellIndex, edge);
                    int otherCellColor;
                    if (otherC != -1)
                    {
                        Cell otherCell = cells[otherC];
                        if (Math.Abs(otherCell.Color) == Math.Abs(cell.Color))
                            continue;
                        otherCellColor = otherCell.Color;
                    }
                    else
                    {
                        if (Math.Abs(cell.Color) == 1)
                            continue;
                        otherCellColor = 1;
                    }
                    List<int>? cellColorSet1 = otherCellColor != 0 ? cellColorSets[Math.Abs(otherCellColor) - 1] : null;
                    if (cellColorSet1 == null)
                    {
                        cellColorSet1 = new List<int>();
                        cellColorSet1.Add(otherC);
                    }
                    List<int> cellColorSet2 = cellColorSets[Math.Abs(cell.Color) - 1];
                    List<int> smallSet;
                    int otherColor;
                    int smallColor;
                    int otherCellIndex;
                    if (cellColorSet1.Count < cellColorSet2.Count)
                    {
                        smallSet = cellColorSet1;
                        smallColor = otherCellColor;
                        otherColor = cell.Color;
                        otherCellIndex = cellIndex;
                    }
                    else
                    {
                        smallSet = cellColorSet2;
                        smallColor = cell.Color;
                        otherColor = otherCellColor;
                        otherCellIndex = otherC;
                    }
                    for (int j = 0; j < smallSet.Count; j++)
                    {
                        int cell3 = smallSet[j];
                        Cell c3 = cells[cell3];
                        for (int i = 0; i < c3.Edges.Count; i++)
                        {
                            int edge2 = c3.Edges[i];
                            if (edge2 == edge)
                                continue;
                            int cell4 = GetAdjacentCell(cell3, edge2);
                            int c4Color;
                            if (cell4 != -1)
                            {
                                Cell c4 = cells[cell4];
                                c4Color = c4.Color;
                            }
                            else
                                c4Color = 1;
                            if (cell4 == otherCellIndex || (c4Color != 0 && Math.Abs(c4Color) == Math.Abs(otherColor)))
                            {
                                // Bingo!
                                int dist = edgeDistances[edge, edge2];
                                if (curDepth + dist >= moves.Length)
                                    continue;
                                if (!AddColorJoinAction(edge, edge2, (c4Color != otherColor) ^ (c3.Color == smallColor), moves, curDepth + dist, edgePairsSeen))
                                    return false;
                            }
                        }
                    }
                }
            }
            return true;
        }

        private bool GatherCellColoringEdgeColoringMovesForEdgeColorChange(Edge edge, List<IAction>[] moves, int curDepth, int edgeIndex)
        {
            if (useCellColoring && useColoring)
            {
                int[,] edgeDistances = this.edgeDistances;
                int cell1 = -1;
                int cell1Color = 0;
                int cell2 = -1;
                int cell2Color = 1;
                for (int i2 = 0; i2 < edge.Cells.Count; i2++)
                {
                    if (cell1 == -1)
                    {
                        cell1 = edge.Cells[i2];
                        cell1Color = cells[cell1].Color;
                    }
                    else
                    {
                        cell2 = edge.Cells[i2];
                        cell2Color = cells[cell2].Color;
                    }
                }
                // If they are the same color we're going to get the inner edge in a second, and that will flow everything as the edges are set.
                // Therefore no point doing this very expensive operation.
                if (cell1Color == cell2Color || cell1Color == -cell2Color)
                    return true;
                List<int> colorSet = colorSets[Math.Abs(edge.Color) - 1];
                for (int i = 0; i < colorSet.Count; i++)
                {
                    int edgeIndex2 = colorSet[i];
                    if (edgeIndex2 == edgeIndex)
                        continue;
                    int dist = edgeDistances[edgeIndex, edgeIndex2];
                    if (dist + curDepth >= moves.Length)
                        continue;
                    Edge edge2 = edges[edgeIndex2];
                    int cell3 = -1;
                    int cell3Color = 0;
                    int cell4 = -1;
                    int cell4Color = 1;
                    List<int> edge2Cells = edge2.Cells;
                    int jMax = edge2Cells.Count;
                    for (int j = 0; j < jMax; j++)
                    {
                        if (cell3 == -1)
                        {
                            cell3 = edge2Cells[j];
                            cell3Color = cells[cell3].Color;
                        }
                        else
                        {
                            cell4 = edge2Cells[j];
                            cell4Color = cells[cell4].Color;
                        }
                    }
                    // we have our 4 edges, what can we do with them... 2x2 we'll try it all...
                    if (cell1 == cell3 || (cell1Color != 0 && (cell1Color == cell3Color || cell1Color == -cell3Color)))
                    {
                        // Real cell overlap.
                        if (cell2 != cell4 && (cell2Color == 0 || (cell2Color !=cell4Color && cell2Color != -cell4Color)))
                        {
                            if (cell2 == -1)
                            {
                                if (!AddCellColorJoinAction(cell4, cell2, (edge.Color != edge2.Color) ^ (cell1Color == cell3Color), moves, curDepth + dist, cellPairsSeen))
                                    return false;
                            }
                            else
                            {
                                if (!AddCellColorJoinAction(cell2, cell4, (edge.Color != edge2.Color) ^ (cell1Color == cell3Color), moves, curDepth + dist, cellPairsSeen))
                                    return false;
                            }
                        }
                        else if ((edge.Color != edge2.Color) ^ (cell2Color == cell4Color) ^ (cell1Color == cell3Color))
                            return false;
                    }
                    if (cell1 == cell4 || (cell1Color != 0 && (cell1Color == cell4Color || cell1Color == -cell4Color)))
                    {
                        // Real cell overlap.
                        if (cell2 != cell3 && (cell2Color == 0 || (cell2Color != cell3Color && cell2Color != -cell3Color)))
                        {
                            if (!AddCellColorJoinAction(cell3, cell2, (edge.Color != edge2.Color) ^ (cell1Color == cell4Color), moves, curDepth + dist, cellPairsSeen))
                                return false;
                        }
                        else if ((edge.Color != edge2.Color) ^ (cell2Color == cell3Color) ^ (cell1Color == cell4Color))
                            return false;
                    }
                    if (cell2 == cell3 || (cell2Color != 0 && (cell2Color == cell3Color || cell2Color == -cell3Color)))
                    {
                        // Real cell overlap.
                        if (cell1 != cell4 && (cell1Color == 0 || (cell1Color != cell4Color && cell1Color != -cell4Color)))
                        {
                            if (!AddCellColorJoinAction(cell1, cell4, (edge.Color != edge2.Color) ^ (cell2Color == cell3Color), moves, curDepth + dist, cellPairsSeen))
                                return false;
                        }
                        else if ((edge.Color != edge2.Color) ^ (cell1Color == cell4Color) ^ (cell2Color == cell3Color))
                            return false;
                    }
                    if (cell2 == cell4 || (cell2Color != 0 && (cell2Color == cell4Color || cell2Color == -cell4Color)))
                    {
                        // Real cell overlap.
                        if (cell1 != cell3 && (cell1Color == 0 || (cell1Color != cell3Color && cell1Color != -cell3Color)))
                        {
                            if (!AddCellColorJoinAction(cell1, cell3, (edge.Color != edge2.Color) ^ (cell2Color == cell4Color), moves, curDepth + dist, cellPairsSeen))
                                return false;
                        }
                        else if ((edge.Color != edge2.Color) ^ (cell1Color == cell3Color) ^ (cell2Color == cell4Color))
                            return false;
                    }
                }
            }
            return true;
        }

        private bool JoinAllLooseCellColors(List<IAction>[] moves, int curDepth, List<int> otherCells, int color, bool pos, int cSpec)
        {
            foreach (int c2 in otherCells)
            {
                if (c2 != -1)
                {
                    Cell c = cells[c2];
                    if (c.Color != color && c.Color != -color)
                    {
                        if (!AddCellColorJoinAction(c2, cSpec, pos, moves, curDepth + 1, cellPairsSeen))
                            return false;
                    }
                }
                else
                {
                    if (color != 1 && -color != -1)
                    {
                        if (!AddCellColorJoinAction(cSpec, c2, pos, moves, curDepth + 1, cellPairsSeen))
                            return false;
                    }
                }
            }
            return true;
        }

        private void FindCellOfSpecificColor(List<int> otherCells, int targetColor, out bool pos, out int c2)
        {
            pos = false;
            c2 = -2;
            foreach (int cell in otherCells)
            {
                if (cell == -1)
                {
                    if (targetColor == 1 || targetColor == -1)
                    {
                        c2 = cell;
                        pos = targetColor == 1;
                        return;
                    }
                }
                else
                {
                    Cell c = cells[cell];
                    if (c.Color == targetColor || c.Color == -targetColor)
                    {
                        c2 = cell;
                        pos = c.Color == targetColor;
                        return;
                    }
                }
            }
        }

        private bool GatherInteractForcedMoves(Intersection inters, List<IAction>[] moves, int curDepth, EdgeState[] edgesSeen, TriState[,] edgePairsSeen, EdgePairRestriction[,] edgeRestrictsSeen)
        {
            if (UseColoring || UseEdgeRestricts)
            {
                // setup coloring and blast away.
                int[] numbering = new int[inters.Edges.Count];
                int[] edgeNumber = new int[inters.Edges.Count];
                int maxExistColor = 0;
                for (int i = 0; i < edgeNumber.Length; i++)
                {
                    edgeNumber[i] = inters.Edges[i];
                    if (UseColoring)
                    {
                        Edge e = edges[edgeNumber[i]];
                        if (e.Color != 0)
                            numbering[i] = e.Color;
                        int colorNum = Math.Abs(e.Color);
                        if (colorNum > maxExistColor)
                            maxExistColor = colorNum;
                    }
                }
                for (int i = 0; i < numbering.Length; i++)
                    if (numbering[i] == 0)
                        numbering[i] = i + 1 + maxExistColor;

                // We're colored, ready to rumble.
                List<int> target = new List<int>();
                if (inters.Type == IntersType.Unknown)
                {
                    target.Add(0);
                }
                target.Add(2);
                // TODO: add the implied edge restricttions/coloring if the inters.Type is not unknown.  AKA replace this with an Advanced Options version.
                if (!GatherFromColoringOptions(moves, curDepth, edgesSeen, numbering, edgeNumber, target, edgePairsSeen, edgeRestrictsSeen))
                    return false;
            }
            return true;
        }


        List<int[]> smallNumberingArrays = new List<int[]>();
        List<int[]> smallEdgeNumberArrays = new List<int[]>();
        List<int> intersTargets = new List<int> { 0, 2 };
        List<List<int>> smallSingleValueLists = new List<List<int>>();

        private bool GatherCellPairForcedMoves(Edge edge, List<IAction>[] moves, int curDepth, EdgeState[] edgesSeen, TriState[,] edgePairsSeen, EdgePairRestriction[,] edgeRestrictsSeen)
        {
            // TODO: replace and rename to intersection pairs.
            /*
            if (UseCellPairs || (topLevel && UseCellPairsTopLevel))
            {
                Cell cell1 = cells[edge.Cells[0]];
                Cell cell2 = null;
                if (edge.Cells.Count > 1)
                    cell2 = cells[edge.Cells[1]];
                Intersection inters1 = intersections[edge.Intersections[0]];
                Intersection inters2 = intersections[edge.Intersections[1]];
                int edgeCount = cell1.Edges.Count;
                // This assumes flat 2 dimensions.
                if (cell2 != null)
                {
                    edgeCount += cell2.Edges.Count - 1;
                    edgeCount += inters1.Edges.Count - 3;
                    edgeCount += inters2.Edges.Count - 3;
                }
                else
                {
                    edgeCount += inters1.Edges.Count - 2;
                    edgeCount += inters2.Edges.Count - 2;
                }
                if (edgeCount >= smallEdgeNumberArrays.Count)
                {
                    for (int i = smallEdgeNumberArrays.Count; i <= edgeCount; i++)
                    {
                        smallNumberingArrays.Add(new int[i]);
                        smallEdgeNumberArrays.Add(new int[i]);
                    }
                }
                int[] numbering = smallNumberingArrays[edgeCount];
                Array.Clear(numbering, 0, numbering.Length);
                int[] edgeNumber = smallEdgeNumberArrays[edgeCount];
                int maxExistColor = 0;
                int counter = 0;
                List<KeyValuePair<int, List<int>>> targets = new List<KeyValuePair<int, List<int>>>();
                int inters1Indexes = 0;
                foreach (int edgeIndex in inters1.Edges)
                {
                    edgeNumber[counter] = edgeIndex;
                    if (UseColoring)
                    {
                        Edge e = edges[edgeNumber[counter]];
                        if (e.Color != 0)
                            numbering[counter] = e.Color;
                        int colorNum = Math.Abs(e.Color);
                        if (colorNum > maxExistColor)
                            maxExistColor = colorNum;
                    }
                    inters1Indexes |= 1 << counter;
                    counter++;
                }
                targets.Add(new KeyValuePair<int, List<int>>(inters1Indexes, intersTargets));
                int inters2Indexes = 0;
                foreach (int edgeIndex in inters2.Edges)
                {
                    Edge e2 = edges[edgeIndex];
                    if (e2 == edge)
                    {
                        inters2Indexes |= 1 << Array.IndexOf(edgeNumber, edgeIndex, 0, counter);
                        continue;
                    }
                    edgeNumber[counter] = edgeIndex;
                    if (UseColoring)
                    {
                        Edge e = edges[edgeNumber[counter]];
                        if (e.Color != 0)
                            numbering[counter] = e.Color;
                        int colorNum = Math.Abs(e.Color);
                        if (colorNum > maxExistColor)
                            maxExistColor = colorNum;
                    }
                    inters2Indexes |= 1 << counter;
                    counter++;
                }
                targets.Add(new KeyValuePair<int, List<int>>(inters2Indexes, intersTargets));
                int cell1Indexes = 0;
                foreach (int edgeIndex in cell1.Edges)
                {
                    Edge e2 = edges[edgeIndex];
                    if (e2 == edge || e2.Intersections[0] == edge.Intersections[0] ||
                        e2.Intersections[0] == edge.Intersections[1] ||
                        e2.Intersections[1] == edge.Intersections[0] ||
                        e2.Intersections[1] == edge.Intersections[1])
                    {
                        cell1Indexes |= 1 << Array.IndexOf(edgeNumber, edgeIndex, 0, counter);
                        continue;
                    }
                    edgeNumber[counter] = edgeIndex;
                    if (UseColoring)
                    {
                        Edge e = edges[edgeNumber[counter]];
                        if (e.Color != 0)
                            numbering[counter] = e.Color;
                        int colorNum = Math.Abs(e.Color);
                        if (colorNum > maxExistColor)
                            maxExistColor = colorNum;
                    }
                    cell1Indexes |= 1 << counter;
                    counter++;
                }
                if (cell1.TargetCount >= 0)
                {
                    if (cell1.TargetCount >= smallSingleValueLists.Count)
                    {
                        for (int i = smallSingleValueLists.Count; i <= cell1.TargetCount; i++)
                        {
                            smallSingleValueLists.Add(new List<int> { i });
                        }
                    }
                    targets.Add(new KeyValuePair<int, List<int>>(cell1Indexes, smallSingleValueLists[cell1.TargetCount]));
                }
                if (cell2 != null)
                {
                    int cell2Indexes = 0;
                    foreach (int edgeIndex in cell2.Edges)
                    {
                        Edge e2 = edges[edgeIndex];
                        if (e2 == edge || e2.Intersections[0] == edge.Intersections[0] ||
                            e2.Intersections[0] == edge.Intersections[1] ||
                            e2.Intersections[1] == edge.Intersections[0] ||
                            e2.Intersections[1] == edge.Intersections[1])
                        {
                            cell2Indexes |= 1 << Array.IndexOf(edgeNumber, edgeIndex, 0, counter);
                            continue;
                        }
                        edgeNumber[counter] = edgeIndex;
                        if (UseColoring)
                        {
                            Edge e = edges[edgeNumber[counter]];
                            if (e.Color != 0)
                                numbering[counter] = e.Color;
                            int colorNum = Math.Abs(e.Color);
                            if (colorNum > maxExistColor)
                                maxExistColor = colorNum;
                        }
                        cell2Indexes |= 1 << counter;
                        counter++;
                    }
                    if (cell2.TargetCount >= 0)
                    {
                        if (cell2.TargetCount >= smallSingleValueLists.Count)
                        {
                            for (int i = smallSingleValueLists.Count; i <= cell2.TargetCount; i++)
                            {
                                smallSingleValueLists.Add(new List<int> { i });
                            }
                        }
                        targets.Add(new KeyValuePair<int, List<int>>(cell2Indexes, smallSingleValueLists[cell2.TargetCount]));
                    }
                }
                for (int i = 0; i < numbering.Length; i++)
                    if (numbering[i] == 0)
                        numbering[i] = i + 1 + maxExistColor;
                if (!GatherFromAdvancedOptions(moves, curDepth, edgesSeen, numbering, edgeNumber, targets, edgePairsSeen, edgeRestrictsSeen))
                    return false;
            }
             */
            return true;
        }



        List<uint> edgeRestrictPattern = new List<uint>();
        List<uint> edgeRestrictMask = new List<uint>();

        private bool GatherFromColoringOptions(List<IAction>[] moves, int curDepth, EdgeState[] edgesSeen, int[] numberingFull, int[] edgeNumber, List<int> targets_in, TriState[,] edgePairsSeen, EdgePairRestriction[,] edgeRestrictsSeen)
        {
            List<KeyValuePair<uint, List<int>>> targets = new List<KeyValuePair<uint, List<int>>> { new KeyValuePair<uint, List<int>>((1u << numberingFull.Length) - 1, targets_in) };
            return GatherFromAdvancedOptions(moves, curDepth, edgesSeen, numberingFull, edgeNumber, targets, edgePairsSeen, edgeRestrictsSeen);
        }


        private struct SuccessLookup
        {
            public uint[] success;
            public int curNumber;

            public override int GetHashCode()
            {
                int hashcode = curNumber;
                for (int i = 0; i < success.Length; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= (int)success[i];
                }
                return hashcode;
            }

            public override bool Equals(object? obj)
            {
                if (obj is not SuccessLookup other)
                    return false;
                if (other.curNumber != curNumber)
                    return false;
                if (other.success.Length != success.Length)
                    return false;
                for (int i = 0; i < success.Length; i++)
                {
                    if (other.success[i] != success[i])
                        return false;
                }
                return true;
              
            }
        }

        Dictionary<SuccessLookup, int[,]> successLookup = new Dictionary<SuccessLookup, int[,]>();

        private int[,] GetMaps(uint[] success, int curNumber)
        {
            SuccessLookup key = new SuccessLookup{success=success, curNumber=curNumber};
            if (!successLookup.TryGetValue(key, out int[,]? map) || map == null)
            {
                map = new int[curNumber, curNumber];
                for (int i = 0; i < curNumber - 1; i++)
                {
                    uint maskI = 1u << i;
                    for (int j = i + 1; j < curNumber; j++)
                    {
                        uint maskJ = 1u << j;
                        int result = 15;
                        for (int k = 0; k < success.Length; k++)
                        {
                            uint val = success[k];
                            if ((val & maskI) != 0 && (val & maskJ) != 0)
                                result &= 14;
                            if ((val & maskI) != 0 && (val & maskJ) == 0)
                                result &= 13;
                            if ((val & maskI) == 0 && (val & maskJ) != 0)
                                result &= 11;
                            if ((val & maskI) == 0 && (val & maskJ) == 0)
                                result &= 7;
                        }
                        map[i, j] = result;
                    }
                }
                successLookup[key] = map;
            }
            return map;
        }

        private bool GatherFromAdvancedOptions(List<IAction>[] moves, int curDepth, EdgeState[] edgesSeen, int[] numberingFull, int[] edgeNumber, List<KeyValuePair<uint, List<int>>> targets, TriState[,] edgePairsSeen, EdgePairRestriction[,] edgeRestrictsSeen)
        {
            int[] baseLine = new int[targets.Count];
            int[] numberingCleared = new int[numberingFull.Length];
            Array.Copy(numberingFull, numberingCleared, numberingFull.Length);
            for (int i = 0; i < edgeNumber.Length; i++)
            {
                Edge e = edges[edgeNumber[i]];
                if (e.State == EdgeState.Empty)
                    continue;
                int index = numberingCleared[i];
                if (e.State != EdgeState.Filled)
                {
                    index = -index;
                }
                if (index != 0)
                {
                    for (int j = 0; j < numberingCleared.Length; j++)
                    {
                        if (numberingCleared[j] == index)
                        {
                            for (int k = 0; k < targets.Count; k++)
                            {
                                if ((targets[k].Key & (1u << j)) != 0)
                                    baseLine[k]++;
                            }
                            numberingCleared[j] = 0;
                        }
                        if (numberingCleared[j] == -index)
                            numberingCleared[j] = 0;
                    }
                }
            }
            int curNumber = 1;
            int[] numbering = new int[numberingCleared.Length];
            for (int i = 0; i < numberingCleared.Length; i++)
            {
                if (numberingCleared[i] == 0 || numbering[i] != 0)
                    continue;
                int cur = numberingCleared[i];
                for (int j = 0; j < numberingCleared.Length; j++)
                {
                    if (numberingCleared[j] == cur)
                        numbering[j] = curNumber;
                    else if (numberingCleared[j] == -cur)
                        numbering[j] = -curNumber;
                }
                curNumber++;
            }
            if (curNumber == 1)
                return true;
            if (useEdgeRestricts)
            {
                MapEdgeRestrictions(edgeNumber, numberingCleared, numbering);
            }
            // make cur number the number of unique numbers for use below rather than making a new variable :P
            curNumber--;
            List<int[]>? result = RetrieveActions(targets, baseLine, curNumber, numbering);
            if (result == null)
                return false;
            return ProcessRetrievedActions(moves, curDepth, edgesSeen, edgeNumber, edgePairsSeen, edgeRestrictsSeen, result);
      }

        private bool ProcessRetrievedActions(List<IAction>[] moves, int curDepth, EdgeState[] edgesSeen, int[] edgeNumber, TriState[,] edgePairsSeen, EdgePairRestriction[,] edgeRestrictsSeen, List<int[]> result)
        {
            for (var i = 0; i < result.Count; i++)
            {
                int[] action = result[i];
                if (action[0] == 0)
                {
                    int j = action[1];
                    bool reallyFilled = action[2] == 1;
                    Edge e = edges[edgeNumber[j]];
                    if (e.State == EdgeState.Empty)
                    {
                        if (!AddSetAction(edgeNumber[j], (reallyFilled ? EdgeState.Filled : EdgeState.Excluded), moves,
                            curDepth + 1, edgesSeen))
                            return false;
                    }
                }
                else if (action[0] == 1)
                {
                    int n = action[1];
                    int m = action[2];
                    bool combined = action[3] == 1;
                    if (!AddColorJoinAction(edgeNumber[n], edgeNumber[m], combined, moves, curDepth + 1, edgePairsSeen))
                        return false;
                }
                else if (action[0] == 2)
                {
                    int n = action[1];
                    int m = action[2];
                    EdgePairRestriction restr = action[3] == 1
                        ? EdgePairRestriction.NotBoth
                        : EdgePairRestriction.NotNeither;
                    if (!AddEdgeRestrictAction(edgeNumber[n], edgeNumber[m], restr, moves, curDepth + 1,
                        edgeRestrictsSeen, edgePairsSeen))
                        return false;
                }
            }
            return true;
        }

        private struct PatternLookup
        {
            public List<KeyValuePair<uint, List<int>>> targets; 
            public int[] baseLine; 
            public int curNumber; 
            public int[] numbering;
            public List<uint> edgePatterns;
            public List<uint> edgeMasks;
            public override int GetHashCode()
            {
                int hashcode = 0;
                for (int i = 0; i < targets.Count; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^=  (int)targets[i].Key;
                    for (int j = 0; j < targets[i].Value.Count; j++)
                    {
                        hashcode += hashcode << 5;
                        hashcode ^= targets[i].Value[j];
                    }
                }
                for (int i = 0; i < baseLine.Length; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= baseLine[i];
                }
                for (int i = 0; i < numbering.Length; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= numbering[i];
                }
                for (int i = 0; i < edgePatterns.Count; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= (int)edgePatterns[i];
                }
                for (int i = 0; i < edgeMasks.Count; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= (int)edgeMasks[i];
                }
                hashcode += hashcode << 5;
                hashcode ^= curNumber;
                return hashcode;
            }
            public override bool Equals(object? obj)
            {
                if (obj is not PatternLookup other)
                    return false;
                if (other.curNumber != curNumber)
                    return false;
                if (other.targets.Count != targets.Count)
                    return false;
                for (int i = 0; i < targets.Count; i++)
                {
                    KeyValuePair<uint, List<int>> entry = targets[i];
                    KeyValuePair<uint, List<int>> otherEntry = other.targets[i];

                    if (otherEntry.Key != entry.Key)
                        return false;
                    if (otherEntry.Value.Count != entry.Value.Count)
                        return false;
                    for (int j = 0; j < entry.Value.Count; j++)
                    {
                        if (otherEntry.Value[j] != entry.Value[j])
                            return false;
                    }
                }
                for (int i = 0; i < baseLine.Length; i++)
                {
                    if(baseLine[i] != other.baseLine[i])
                        return false;
                }
                for (int i = 0; i < numbering.Length; i++)
                {
                    if (numbering[i] != other.numbering[i])
                        return false;
                }
                if (edgePatterns.Count != other.edgePatterns.Count)
                    return false;
                for (int i = 0; i < edgePatterns.Count; i++)
                {
                    if (edgePatterns[i] != other.edgePatterns[i])
                        return false;
                }
                for (int i = 0; i < edgeMasks.Count; i++)
                {
                    if (edgeMasks[i] != other.edgeMasks[i])
                        return false;
                }
                return true;
            }
        }

        Dictionary<PatternLookup, List<int[]>?> patternLookup = new Dictionary<PatternLookup, List<int[]>?>();
        List<uint> emptyList = new List<uint>();

        private List<int[]>? RetrieveActions(List<KeyValuePair<uint, List<int>>> targets, int[] baseLine, int curNumber, int[] numbering)
        {
            List<int[]>? result;
            PatternLookup key;
            if (useEdgeRestricts)
            {
                key = new PatternLookup() { targets = targets, baseLine = baseLine, curNumber = curNumber, numbering = numbering, edgePatterns = edgeRestrictPattern, edgeMasks = edgeRestrictMask };
            }
            else
            {
                key = new PatternLookup() { targets = targets, baseLine = baseLine, curNumber = curNumber, numbering = numbering, edgePatterns = emptyList, edgeMasks = emptyList };
            }
            if (!patternLookup.TryGetValue(key, out result))
            {
                if (useEdgeRestricts)
                {
                    // edgeRestrictPattern and mask are shared data, we must clone them before we store the key.
                    // We clone inside the try get value to avoid the cost of cloning them for lookup.
                    key.edgePatterns = new List<uint>(edgeRestrictPattern);
                    key.edgeMasks = new List<uint>(edgeRestrictMask);
                }

                int[] countsFor = new int[targets.Count * curNumber];
                int[] countsAgainst = new int[targets.Count * curNumber];
                for (int j = 0; j < targets.Count; j++)
                {
                    uint checkPattern = targets[j].Key;
                    for (int i = 0; i < numbering.Length; i++)
                    {
                        if ((checkPattern & (1u << i)) != 0u)
                        {
                            if (numbering[i] > 0)
                                countsFor[j * curNumber + numbering[i] - 1]++;
                            else if (numbering[i] < 0)
                                countsAgainst[j * curNumber - numbering[i] - 1]++;
                        }
                    }
                }
                List<uint> successIn = new List<uint>();

                uint max = 1u << curNumber;
                int[][] targetCounts = new int[targets.Count][];
                for (int j = 0; j < targets.Count; j++)
                {
                    targetCounts[j] = targets[j].Value.ToArray();
                }
                for (uint i = 0; i < max; i++)
                {
                    bool fail = false;
                    for (int k = 0; k < targets.Count; k++)
                    {
                        int total = baseLine[k];
                        for (int j = 0; j < curNumber; j++)
                        {
                            if ((i & (1u << j)) != 0)
                            {
                                total += countsFor[k * curNumber + j];
                            }
                            else
                            {
                                total += countsAgainst[k * curNumber + j];
                            }
                        }
                        int[] targetset = targetCounts[k];
                        int targetCount = targetset.Length;
                        if (targetCount == 1 && targetset[0] != total)
                        {
                            fail = true;
                            break;
                        }
                        else if (targetCount == 2 && targetset[0] != total && targetset[1] != total)
                        {
                            fail = true;
                            break;
                        }
                    }
                    if (!fail)
                        successIn.Add(i);
                }
                if (useEdgeRestricts && edgeRestrictPattern.Count != 0)
                {
                    List<uint> success2 = new List<uint>();
                    for (int j = 0; j < successIn.Count; j++)
                    {
                        uint i = successIn[j];
                        bool fail = false;
                        for (int k = 0; k < edgeRestrictPattern.Count; k++)
                        {
                            if ((i & edgeRestrictMask[k]) == edgeRestrictPattern[k])
                            {
                                fail = true;
                                break;
                            }
                        }
                        if (fail)
                            continue;
                        success2.Add(i);
                    }
                    successIn = success2;
                }
                if (successIn.Count == 0)
                    result = null;
                else
                {
                    result = new List<int[]>();
                    uint[] success = successIn.ToArray();
                    uint set = uint.MaxValue;
                    for (int i = 0; i < success.Length; i++)
                        set &= success[i];
                    uint unset = uint.MaxValue;
                    for (int i = 0; i < success.Length; i++)
                        unset &= ~success[i];
                    for (int i = 0; i < curNumber; i++)
                    {
                        bool filled;
                        if ((set & (1u << i)) != 0)
                        {
                            filled = true;
                        }
                        else if ((unset & (1u << i)) != 0)
                        {
                            filled = false;
                        }
                        else
                            continue;
                        for (int j = 0; j < numbering.Length; j++)
                        {
                            bool follow;
                            if (numbering[j] == i + 1)
                            {
                                follow = true;
                            }
                            else if (-numbering[j] == i + 1)
                            {
                                follow = false;
                            }
                            else
                                continue;
                            bool reallyFilled = !follow ^ filled;
                            result.Add(new int[] { 0, j, reallyFilled ? 1 : 0 });
                        }
                    }
                    int[,]? maps = null;
                    if (UseColoring || UseEdgeRestricts)
                        maps = GetMaps(success, curNumber);
                    if (UseColoring)
                    {
                        for (int i = 0; i < curNumber - 1; i++)
                        {
                            for (int j = i + 1; j < curNumber; j++)
                            {
                                int value = maps![i, j];
                                bool same = ((value & 2) != 0) && ((value & 4) != 0);
                                bool opposite = ((value & 1) != 0) && ((value & 8) != 0);
                                if (same || opposite)
                                {
                                    bool invert = opposite == true;
                                    // i and j can be joinedin the current collapsed numbering.
                                    // unnumbered sections need to be joined, 
                                    for (int n = 0; n < numbering.Length; n++)
                                    {
                                        int val = Math.Abs(numbering[n]);
                                        if (val == i + 1 || val == j + 1)
                                        {
                                            bool pos = numbering[n] > 0;
                                            for (int m = 0; m < numbering.Length; m++)
                                            {
                                                if (m == n)
                                                    continue;
                                                int val2 = Math.Abs(numbering[m]);
                                                if (val2 == j + 1 || val2 == i + 1)
                                                {
                                                    bool pos2 = numbering[m] > 0;
                                                    bool combined = pos2 == pos;
                                                    if (val2 != val)
                                                        combined ^= invert;
                                                    result.Add(new int[] { 1, n, m, combined ? 1 : 0 });
                                                }
                                            }
                                        }
                                    }
                                }

                            }
                        }
                    }
                    // TODO: this is going to rediscover the same edge restrictions we have as input.
                    // Do not add them to the output, as it just inflates output for no good reason.
                    if (UseEdgeRestricts)
                    {
                        for (int i = 0; i < curNumber - 1; i++)
                        {
                            for (int j = i + 1; j < curNumber; j++)
                            {
                                int value = maps![i, j];
                                bool not11 = (value & 1) != 0;
                                bool not00 = (value & 8) != 0;
                                bool not10 = (value & 2) != 0;
                                bool not01 = (value & 4) != 0;
                                int countTrue = 0;
                                if (not11)
                                    countTrue++;
                                if (not00)
                                    countTrue++;
                                if (not01)
                                    countTrue++;
                                if (not10)
                                    countTrue++;
                                if (countTrue == 1)
                                {
                                    for (int n = 0; n < numbering.Length; n++)
                                    {
                                        int val = Math.Abs(numbering[n]);
                                        if (val == i + 1)
                                        {
                                            bool pos = numbering[n] > 0;
                                            for (int m = 0; m < numbering.Length; m++)
                                            {
                                                if (m == n)
                                                    continue;
                                                int val2 = Math.Abs(numbering[m]);
                                                if (val2 == j + 1)
                                                {
                                                    bool pos2 = numbering[m] > 0;
                                                    if (not11)
                                                    {
                                                        if (pos && pos2)
                                                        {
                                                            result.Add(new int[] { 2, n, m, 1 });
                                                        }
                                                        else if (!pos && !pos2)
                                                        {
                                                            result.Add(new int[] { 2, n, m, 0 });
                                                        }
                                                    }
                                                    else if (not00)
                                                    {
                                                        if (pos && pos2)
                                                        {
                                                            result.Add(new int[] { 2, n, m, 0 });
                                                        }
                                                        else if (!pos && !pos2)
                                                        {
                                                            result.Add(new int[] { 2, n, m, 1 });
                                                        }
                                                    }
                                                    else if (not10)
                                                    {
                                                        if (pos && !pos2)
                                                        {
                                                            result.Add(new int[] { 2, n, m, 1 });
                                                        }
                                                        else if (!pos && pos2)
                                                        {
                                                            result.Add(new int[] { 2, n, m, 0 });
                                                        }
                                                    }
                                                    else if (not01)
                                                    {
                                                        if (!pos && pos2)
                                                        {
                                                            result.Add(new int[] { 2, n, m, 1 });
                                                        }
                                                        else if (pos && !pos2)
                                                        {
                                                            result.Add(new int[] { 2, n, m, 0 });
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                            }
                        }
                    }
                }
                patternLookup[key] = result;
            }
            return result;
        }

        private void MapEdgeRestrictions(int[] edgeNumber, int[] numberingCleared, int[] numbering)
        {
            edgeRestrictPattern.Clear();
            edgeRestrictMask.Clear();
            if (edgePairRestrictions == null)
                return;
            for (int i = 0; i < numberingCleared.Length - 1; i++)
            {
                int firstEdge = edgeNumber[i];
                if (numbering[i] == 0)
                    continue;
                for (int j = i + 1; j < numberingCleared.Length; j++)
                {
                    int secondEdge = edgeNumber[j];
                    if (numbering[j] == 0)
                        continue;
                    EdgePairRestriction restr = edgePairRestrictions[firstEdge, secondEdge];
                    if (restr == EdgePairRestriction.NotBoth)
                    {
                        if (numbering[i] > 0 && numbering[j] > 0)
                        {
                            edgeRestrictMask.Add((1u << (numbering[i] - 1)) | (1u << (numbering[j] - 1)));
                            edgeRestrictPattern.Add((1u << (numbering[i] - 1)) | (1u << (numbering[j] - 1)));
                        }
                        else if (numbering[i] > 0)
                        {
                            if (numbering[i] != -numbering[j])
                            {
                                edgeRestrictMask.Add((1u << (numbering[i] - 1)) | (1u << (-numbering[j] - 1)));
                                edgeRestrictPattern.Add((1u << (numbering[i] - 1)));
                            }
                        }
                        else if (numbering[j] > 0)
                        {
                            if (numbering[i] != -numbering[j])
                            {
                                edgeRestrictMask.Add((1u << (-numbering[i] - 1)) | (1u << (numbering[j] - 1)));
                                edgeRestrictPattern.Add((1u << (numbering[j] - 1)));
                            }
                        }
                        else
                        {
                            edgeRestrictMask.Add((1u << (-numbering[i] - 1)) | (1u << (-numbering[j] - 1)));
                            edgeRestrictPattern.Add(0);
                        }
                    }
                    else if (restr == EdgePairRestriction.NotNeither)
                    {
                        if (numbering[i] > 0 && numbering[j] > 0)
                        {
                            edgeRestrictMask.Add((1u << (numbering[i] - 1)) | (1u << (numbering[j] - 1)));
                            edgeRestrictPattern.Add(0);
                        }
                        else if (numbering[i] > 0)
                        {
                            if (numbering[i] != -numbering[j])
                            {
                                edgeRestrictMask.Add((1u << (numbering[i] - 1)) | (1u << (-numbering[j] - 1)));
                                edgeRestrictPattern.Add((1u << (-numbering[j] - 1)));
                            }
                        }
                        else if (numbering[j] > 0)
                        {
                            if (numbering[i] != -numbering[j])
                            {
                                edgeRestrictMask.Add((1u << (-numbering[i] - 1)) | (1u << (numbering[j] - 1)));
                                edgeRestrictPattern.Add((1u << (-numbering[i] - 1)));
                            }
                        }
                        else
                        {
                            edgeRestrictMask.Add((1u << (-numbering[i] - 1)) | (1u << (-numbering[j] - 1)));
                            edgeRestrictPattern.Add((1u << (-numbering[i] - 1)) | (1u << (-numbering[j] - 1)));
                        }
                    }
                }
            }
        }



    }
}