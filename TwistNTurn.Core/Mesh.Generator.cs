using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;

namespace TwistNTurn
{
    public partial class Mesh
    {

        public bool CreatingLoop = false;

        private bool TryInitialLoopAttempt(Random rnd, List<IAction> backup, int targetCount,
            ref long countSum, ref int tries, ref int loopTries)
        {
            FullClear();
            backup.Clear();
            int start = rnd.Next(intersections.Count);
            Intersection inters = intersections[start];
            int edgeIntersIndex = rnd.Next(inters.Edges.Count);
            int edgeIndex = inters.Edges[edgeIntersIndex];
            Perform(edgeIndex, EdgeState.Filled, backup, 0);
            _ = CreateLoop(rnd, start, start, edgeIntersIndex);
            int count = 0;
            for (var index = 0; index < edges.Count; index++)
            {
                if (edges[index].State == EdgeState.Filled)
                    count++;
            }
            countSum += count;
            tries++;
            if (count < targetCount || RateBoringness() > GenerateBoringFraction)
            {
                if (loopTries > 100 && count >= countSum / tries)
                {
                    if (TryExpandLoop(rnd, targetCount - count))
                        return true;
                }
                loopTries++;
                return loopTries > 1000;
            }
            return true;
        }

        private void GenerateInitialLoop(Random rnd, List<IAction> backup)
        {
            int targetCount = (int)Math.Floor(intersections.Count * generateLengthFraction);
            long countSum = 0;
            int tries = 0;
            int loopTries = 0;
            while (!TryInitialLoopAttempt(rnd, backup, targetCount, ref countSum, ref tries, ref loopTries))
            {
            }
        }

        private async Task GenerateInitialLoopAsync(Random rnd, List<IAction> backup, CancellationToken cancellationToken = default)
        {
            int targetCount = (int)Math.Floor(intersections.Count * generateLengthFraction);
            long countSum = 0;
            int tries = 0;
            int loopTries = 0;
            long lastYieldTicks = Stopwatch.GetTimestamp();

            while (!AbortPrune && !cancellationToken.IsCancellationRequested)
            {
                long currentTicks = Stopwatch.GetTimestamp();
                double elapsedMs = (currentTicks - lastYieldTicks) * 1000.0 / Stopwatch.Frequency;
                if (elapsedMs >= 16.0)
                {
                    lastYieldTicks = currentTicks;
                    await Task.Delay(1, cancellationToken);
                    if (AbortPrune || cancellationToken.IsCancellationRequested)
                        return;
                }

                if (TryInitialLoopAttempt(rnd, backup, targetCount, ref countSum, ref tries, ref loopTries))
                    break;
            }
        }

        private bool StepLoopOptimization(int iteration)
        {
            UpdateTwists();
            if (!MakeLoopBetter())
                return false;
            ClearTwists();
            if (iteration >= 7)
                UpdateTwists();
            return true;
        }

        public void Generate()
        {
            AbortPrune = false;
            bool done = false;
            List<IAction> backup = [];
            Random rnd = Random.Shared;
            while (!done)
            {
                try
                {
                    CreatingLoop = true;
                    GenerateInitialLoop(rnd, backup);
                    for (int i = 0; i < 8; i++)
                    {
                        if (!StepLoopOptimization(i))
                            break;
                    }
                }
                finally
                {
                    CreatingLoop = false;
                }
                Clear();
                done = TryPruneCounts();
            }
        }

        private bool MakeLoopBetter()
        {
            int filledEdge = -1;
            int filledCount = 0;
            for (int i = 0; i < edges.Count; i++)
            {
                if (edges[i].State == EdgeState.Filled)
                {
                    filledEdge = i;
                    filledCount++;
                }
            }
            List<int> edgesSeen = [with(filledCount + 1)];
            List<int> intersSeen = [with(filledCount)];
            List<IntersType> intersTypes = [with(filledCount)];
            List<bool> turned = [with(filledCount)];
            int nextInter = edges[filledEdge].Intersections[0];
            edgesSeen.Add(filledEdge);
            for (int i = 0; i < filledCount; i++)
            {
                intersSeen.Add(nextInter);
                Intersection nextIntersection = intersections[nextInter];
                intersTypes.Add(nextIntersection.Type);
                int nextEdge = GetNextEdge(nextInter, filledEdge);
                int oppositeEdge = nextIntersection.Opposite[nextIntersection.Edges.IndexOf(filledEdge)];
                turned.Add(nextEdge != oppositeEdge);
                edgesSeen.Add(nextEdge);
                filledEdge = nextEdge;
                nextInter = GetOtherInters(filledEdge, nextInter);
            }
            edgesSeen.RemoveAt(edgesSeen.Count - 1);
            List<int> candidates = [];
            for (int i = 0; i < filledCount; i++)
            {
                if (intersTypes[i] != IntersType.Unknown)
                    continue;
                else if (intersTypes[(i + 1) % filledCount] == IntersType.Turn)
                    continue;
                else if (intersTypes[(i + filledCount - 1) % filledCount] != IntersType.Unknown)
                    continue;
                else if (intersTypes[(i + filledCount - 2) % filledCount] == IntersType.Turn)
                    continue;
                // if turned, followed by twist, followed by not turned, in either direction, continue;
                else if (turned[i] && intersTypes[(i + 1) % filledCount] == IntersType.Twist && !turned[(i + 2) % filledCount])
                    continue;
                else if (turned[(i + filledCount - 1) % filledCount] && intersTypes[(i + filledCount - 2) % filledCount] == IntersType.Twist && !turned[(i + filledCount - 3) % filledCount])
                    continue;
                candidates.Add(edgesSeen[i]);
            }
            HashSet<int> candidateCheck = [with(candidates.Count), ..candidates];
            bool changed = false;
            foreach (int edge in candidates)
            {
                changed |= TryPop(edge, candidateCheck);
            }
            return changed;
        }


        private bool TryPop(int edgeToPop, HashSet<int> candidateCheck)
        {
            Edge eToPop = edges[edgeToPop];
            for (var i = 0; i < eToPop.Cells.Count; i++)
            {
                int c = eToPop.Cells[i];
                Cell cell = cells[c];
                if (!IsPopable(cell, candidateCheck))
                    continue;
                int oldCount = cell.FilledCount;
                for (var index = 0; index < cell.Edges.Count; index++)
                {
                    int edge = cell.Edges[index];
                    if (edges[edge].State == EdgeState.Filled)
                    {
                        Perform(new UnsetAction(this, edge), new List<IAction>(), 0, null);
                    }
                    else
                    {
                        if (edges[edge].State == EdgeState.Excluded)
                        {
                            Perform(new UnsetAction(this, edge), new List<IAction>(), 0, null);
                        }
                        Perform(edge, EdgeState.Filled, new List<IAction>(), 0);
                        Edge e = edges[edge];
                        foreach (int interIndex in e.Intersections)
                        {
                            Intersection inter = intersections[interIndex];
                            if (inter.FilledCount == 2)
                            {
                                for (var index1 = 0; index1 < inter.Edges.Count; index1++)
                                {
                                    int edgeIndex = inter.Edges[index1];
                                    Edge otherEdge = edges[edgeIndex];
                                    if (otherEdge.State == EdgeState.Empty)
                                    {
                                        Perform(edgeIndex, EdgeState.Excluded, new List<IAction>(), 0);
                                    }
                                }
                            }
                        }
                    }
                }
                return true;
            }
            return false;
        }

        private bool IsPopable(Cell cell, HashSet<int> candidateCheck)
        {
            if (cell.FilledCount == 0 || cell.FilledCount == cell.Edges.Count)
                return false;
            if (GlancingTouch(cell))
                return false;
            if (!ContiguousFilledEdge(cell))
                return false;
            if (!AllCandidates(cell, candidateCheck))
                return false;
            return true;
        }

        private bool AllCandidates(Cell cell, HashSet<int> candidateCheck)
        {
            for (var index = 0; index < cell.Edges.Count; index++)
            {
                int edge = cell.Edges[index];
                if (edges[edge].State == EdgeState.Filled && !candidateCheck.Contains(edge))
                    return false;
            }
            return true;
        }

        bool[]? boringEdges;

        private double RateBoringness()
        {
            if (boringEdges == null)
                boringEdges = new bool[edges.Count];
            else
                Array.Clear(boringEdges, 0, boringEdges.Length);
            int boringCount = 0;
            int total = edges.Count;
            edgeLoop: for (int index = 0; index < edges.Count; index++)
            {
                Edge e = edges[index];
                if (e.State != EdgeState.Empty)
                    continue;
                foreach (int interIndex in e.Intersections)
                {
                    Intersection inter = intersections[interIndex];
                    for (var i = 0; i < inter.Edges.Count; i++)
                    {
                        if (edges[inter.Edges[i]].State != EdgeState.Empty)
                            continue edgeLoop;
                    }
                }
                boringCount++;
                boringEdges[index] = true;
            }
            return (double)boringCount / (double)total;
        }

        private bool TryExpandLoop(Random rnd, int targetAdditional)
        {
            int lastTarget = 0;
            while (lastTarget != targetAdditional)
            {
                lastTarget = targetAdditional;
                List<int> cellsToExpand = [];
                for (int tries = 0; tries < 2; tries++)
                {
                    for (int i = 0; i < cells.Count; i++)
                    {
                        if (IsExpandable(cells[i]))
                        {
                            if (targetAdditional > 0 || tries != 0 || NearBoringEdge(cells[i]))
                                cellsToExpand.Add(i);
                        }
                    }
                    if (cellsToExpand.Count > 0 || targetAdditional > 0)
                        break;
                }
                rnd.Shuffle(CollectionsMarshal.AsSpan(cellsToExpand));
                int index = 0;
                while ((targetAdditional > 0 || RateBoringness() > GenerateBoringFraction) && index < cellsToExpand.Count)
                {
                    targetAdditional -= Expand(cells[cellsToExpand[index]]);
                    if (rnd.Next(2) == 0)
                        break;
                    index++;
                }
                if (targetAdditional <= 0 && RateBoringness() <= GenerateBoringFraction)
                    return true;
            }
            return false;
        }

        private bool NearBoringEdge(Cell cell)
        {
            if (boringEdges == null)
                RateBoringness();
            for (var i = 0; i < cell.Intersections.Count; i++)
            {
                int interIndex = cell.Intersections[i];
                Intersection inter = intersections[interIndex];
                for (var index = 0; index < inter.Edges.Count; index++)
                {
                    int edge = inter.Edges[index];
                    if (boringEdges![edge])
                        return true;
                }
            }
            return false;
        }

        private bool IsExpandable(Cell cell)
        {
            if (cell.FilledCount == 0)
                return false;
            if (cell.FilledCount >= (cell.Edges.Count + 1) / 2)
                return false;
            if (GlancingTouch(cell))
                return false;
            if (!ContiguousFilledEdge(cell))
                return false;
            return true;
        }

        private bool GlancingTouch(Cell cell)
        {
            for (var index1 = 0; index1 < cell.Intersections.Count; index1++)
            {
                int interIndex = cell.Intersections[index1];
                Intersection inter = intersections[interIndex];
                if (inter.FilledCount != 2)
                    continue;
                bool found = false;
                for (var i = 0; i < inter.Edges.Count; i++)
                {
                    int edgeIndex = inter.Edges[i];
                    Edge edge = edges[edgeIndex];
                    if (edge.State == EdgeState.Filled)
                    {
                        for (var index = 0; index < edge.Cells.Count; index++)
                        {
                            int cellIndex = edge.Cells[index];
                            if (cells[cellIndex] == cell)
                                found = true;
                        }
                    }
                }
                if (!found)
                    return true;
            }
            return false;
        }

        private bool ContiguousFilledEdge(Cell cell)
        {
            if (cell.FilledCount < 2 || cell.FilledCount >= cell.Edges.Count-1)
                return true;
            bool[] filled = new bool[cell.Edges.Count];
            // Relies on intersections being in clock/anticlock order.
            for (int i = 0; i < cell.Intersections.Count; i++)
            {
                int edge = GetEdgeJoining(cell.Intersections[i], cell.Intersections[(i + 1) % cell.Intersections.Count]);
                filled[i] = edges[edge].State == EdgeState.Filled;
            }
            int transitions =0;
            for (int i = 0; i < filled.Length; i++)
            {
                if (filled[i] != filled[(i + 1) % filled.Length])
                {
                    transitions++;
                }
            }
            return transitions <= 2;
        }

        private int Expand(Cell cell)
        {
            if (!IsExpandable(cell))
                return 0;
            int oldCount = cell.FilledCount;
            for (var index = 0; index < cell.Edges.Count; index++)
            {
                int edge = cell.Edges[index];
                if (edges[edge].State == EdgeState.Filled)
                {
                    Perform(new UnsetAction(this, edge), new List<IAction>(), 0, null);
                }
                else
                {
                    if (edges[edge].State == EdgeState.Excluded)
                    {
                        Perform(new UnsetAction(this, edge), new List<IAction>(), 0, null);
                    }
                    Perform(edge, EdgeState.Filled, new List<IAction>(), 0);
                    Edge e = edges[edge];
                    foreach (int interIndex in e.Intersections)
                    {
                        Intersection inter = intersections[interIndex];
                        if (inter.FilledCount == 2)
                        {
                            for (var i = 0; i < inter.Edges.Count; i++)
                            {
                                int edgeIndex = inter.Edges[i];
                                Edge otherEdge = edges[edgeIndex];
                                if (otherEdge.State == EdgeState.Empty)
                                {
                                    Perform(edgeIndex, EdgeState.Excluded, new List<IAction>(), 0);
                                }
                            }
                        }
                    }
                }
            }
            return cell.FilledCount - oldCount;
        }


        private bool CreateLoop(Random rnd, int start, int prev, int prevIntersEdge)
        {
            Intersection prevInters = intersections[prev];
            Edge prevEdge = edges[prevInters.Edges[prevIntersEdge]];
            Intersection? curInters = null;
            int curIntersIndex = -1;
            for (var i = 0; i < prevEdge.Intersections.Length; i++)
            {
                int intersPos = prevEdge.Intersections[i];
                Intersection other = intersections[intersPos];
                if (other != prevInters)
                {
                    curIntersIndex = intersPos;
                    curInters = other;
                    break;
                }
            }
            if (curInters == null)
                return false;
            if (curIntersIndex == start)
                return true;
            if (!CanReach(start, curIntersIndex))
                return false;
            List<int> availEdges = new List<int>();
            for (int i = 0; i < curInters.Edges.Count; i++)
            {
                Edge edge = edges[curInters.Edges[i]];
                if (edge.State == EdgeState.Empty)
                {
                    availEdges.Add(i);
                }
            }
            if (availEdges.Count == 0)
                return false;
            if (availEdges.Count > 1)
            {
                rnd.Shuffle(CollectionsMarshal.AsSpan(availEdges));
            }
            List<IAction> backup = new List<IAction>();
            for (var i2 = 0; i2 < availEdges.Count; i2++)
            {
                int availEdge = availEdges[i2];
                backup.Clear();
                int edgeIndex = curInters.Edges[availEdge];
                Perform(edgeIndex, EdgeState.Filled, backup, 0);

                for (int i = 0; i < curInters.Edges.Count; i++)
                {
                    int otherEdgeIndex = curInters.Edges[i];
                    Edge edge = edges[otherEdgeIndex];
                    if (edge.State == EdgeState.Empty)
                        Perform(otherEdgeIndex, EdgeState.Excluded, backup, 0);
                }
                if (!CreateLoop(rnd, start, curIntersIndex, availEdge))
                {
                    Unperform(backup);
                }
                else
                    return true;
            }
            return false;

        }

        private bool CanReach(int start, int curIntersIndex)
        {
            bool[] reached = new bool[intersections.Count];
            List<int> newReached = new List<int>();
            reached[curIntersIndex] = true;
            newReached.Add(curIntersIndex);
            for (int i = 0; i < newReached.Count; i++)
            {
                int index = newReached[i];
                Intersection inters = intersections[index];
                for (int j=0; j < inters.Edges.Count; j++)
                {
                    int edgeIndex = inters.Edges[j];
                    Edge edge = edges[edgeIndex];
                    if (edge.State == EdgeState.Empty)
                    {
                        for (int k=0; k < edge.Intersections.Length; k++)
                        {
                            int otherInters = edge.Intersections[k];
                            if (otherInters != index)
                            {
                                if (!reached[otherInters])
                                {
                                    if (otherInters == start)
                                        return true;
                                    reached[otherInters] = true;
                                    newReached.Add(otherInters);
                                }
                            }
                        }
                    }
                }
            }
            return reached[start];
        }

        private void UpdateTwists()
        {
            for (int i = 0; i < intersections.Count; i++)
            {
                Intersection inter = intersections[i];
                if (inter.FilledCount == 2)
                {
                    inter.Type = DetermineType(i);
                }
            }
        }

        private void ClearTwists()
        {
            for (int i = 0; i < intersections.Count; i++)
            {
                Intersection inter = intersections[i];
                inter.Type = IntersType.Unknown;
            }
        }

        public IntersType DetermineType(int inter)
        {
            int next1, next2;
            if (IsStraightWithInters(inter, out next1, out next2))
            {
                return (IsTurn(next1) || IsTurn(next2)) ? IntersType.Twist : IntersType.Unknown;
            }
            if (IsTurnWithInters(inter, out next1, out next2))
            {
                return IsStraight(next1) && IsStraight(next2) ? IntersType.Turn : IntersType.Unknown;
            }
            return IntersType.Unknown;
        }

        private IntersType DetermineTypeAssuming(int inter, int extraEdge)
        {
            int next1, next2;
            if (IsStraightWithIntersAssuming(inter, extraEdge, out next1, out next2))
            {
                return (IsTurnAssuming(next1, extraEdge) || IsTurnAssuming(next2, extraEdge)) ? IntersType.Twist : IntersType.Unknown;
            }
            if (IsTurnWithIntersAssuming(inter, extraEdge, out next1, out next2))
            {
                return IsStraightAssuming(next1, extraEdge) && IsStraightAssuming(next2, extraEdge) ? IntersType.Turn : IntersType.Unknown;
            }
            return IntersType.Unknown;
        }
        private bool IsTurnWithInters(int inter, out int next1, out int next2)
        {
            next1 = -1;
            next2 = -1;
            int nextE1, nextE2;
            bool turn = IsTurn(inter, out nextE1, out nextE2);
            if (turn)
            {
                next1 = GetOtherInters(nextE1, inter);
                next2 = GetOtherInters(nextE2, inter);
            }
            return turn;
        }

        private bool IsTurn(int interIndex, out int nextE1, out int nextE2)
        {
            nextE1 = -1;
            nextE2 = -1;
            Intersection inter = intersections[interIndex];
            if (inter.FilledCount != 2)
                return false;
            for (int i = 0; i < inter.Edges.Count; i++)
            {
                if (edges[inter.Edges[i]].State == EdgeState.Filled)
                {
                    if (inter.Opposite[i] >= 0 && edges[inter.Opposite[i]].State == EdgeState.Filled)
                    {
                        return false;
                    }
                    else
                    {
                        nextE1 = inter.Edges[i];
                        for (int j = i+1; j < inter.Edges.Count; j++)
                        {
                            if (edges[inter.Edges[j]].State == EdgeState.Filled)
                            {
                                nextE2 = inter.Edges[j];
                            }
                        }
                        return true;
                    }
                }
            }
            return false;
        }

        private bool IsTurn(int interIndex)
        {
            int a, b;
            return IsTurn(interIndex, out a, out b);
        }

        private bool IsStraightWithInters(int inter, out int next1, out int next2)
        {
            next1 = -1;
            next2 = -1;
            int nextE1, nextE2;
            bool straight = IsStraight(inter, out nextE1, out nextE2);
            if (straight)
            {
                next1 = GetOtherInters(nextE1, inter);
                next2 = GetOtherInters(nextE2, inter);
            }
            return straight;
        }

        private bool IsStraight(int interIndex, out int nextE1, out int nextE2)
        {
            nextE1 = -1;
            nextE2 = -1;
            Intersection inter = intersections[interIndex];
            if (inter.FilledCount != 2)
                return false;
            for (int i = 0; i < inter.Edges.Count; i++)
            {
                if (edges[inter.Edges[i]].State == EdgeState.Filled)
                {
                    if (inter.Opposite[i] >= 0 && edges[inter.Opposite[i]].State == EdgeState.Filled)
                    {
                        nextE1 = inter.Edges[i];
                        nextE2 = inter.Opposite[i];
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
            return false;
        }

        private bool IsStraight(int interIndex)
        {
            int a, b;
            return IsStraight(interIndex, out a, out b);
        }

        private bool IsTurnWithIntersAssuming(int inter, int extraEdge, out int next1, out int next2)
        {
            next1 = -1;
            next2 = -1;
            int nextE1, nextE2;
            bool turn = IsTurnAssuming(inter, extraEdge, out nextE1, out nextE2);
            if (turn)
            {
                next1 = GetOtherInters(nextE1, inter);
                next2 = GetOtherInters(nextE2, inter);
            }
            return turn;
        }

        private bool IsTurnAssuming(int interIndex, int extraEdge, out int nextE1, out int nextE2)
        {
            nextE1 = -1;
            nextE2 = -1;
            Intersection inter = intersections[interIndex];
            if (inter.FilledCount != 2 && !(inter.FilledCount == 1 && inter.Edges.Contains(extraEdge)))
                return false;
            for (int i = 0; i < inter.Edges.Count; i++)
            {
                if (edges[inter.Edges[i]].State == EdgeState.Filled || inter.Edges[i]==extraEdge)
                {
                    if (inter.Opposite[i] >= 0 && (edges[inter.Opposite[i]].State == EdgeState.Filled || inter.Opposite[i] == extraEdge))
                    {
                        return false;
                    }
                    else
                    {
                        nextE1 = inter.Edges[i];
                        for (int j = i + 1; j < inter.Edges.Count; j++)
                        {
                            if (edges[inter.Edges[j]].State == EdgeState.Filled || inter.Edges[j] == extraEdge)
                            {
                                nextE2 = inter.Edges[j];
                            }
                        }
                        return true;
                    }
                }
            }
            return false;
        }

        private bool IsTurnAssuming(int interIndex, int extraEdge)
        {
            int a, b;
            return IsTurnAssuming(interIndex, extraEdge, out a, out b);
        }

        private bool IsStraightWithIntersAssuming(int inter, int extraEdge, out int next1, out int next2)
        {
            next1 = -1;
            next2 = -1;
            int nextE1, nextE2;
            bool straight = IsStraightAssuming(inter, extraEdge, out nextE1, out nextE2);
            if (straight)
            {
                next1 = GetOtherInters(nextE1, inter);
                next2 = GetOtherInters(nextE2, inter);
            }
            return straight;
        }

        private bool IsStraightAssuming(int interIndex, int extraEdge, out int nextE1, out int nextE2)
        {
            nextE1 = -1;
            nextE2 = -1;
            Intersection inter = intersections[interIndex];
            if (inter.FilledCount != 2 && !(inter.FilledCount == 1 && inter.Edges.Contains(extraEdge)))
                return false;
            for (int i = 0; i < inter.Edges.Count; i++)
            {
                if (edges[inter.Edges[i]].State == EdgeState.Filled || inter.Edges[i] == extraEdge)
                {
                    if (inter.Opposite[i] >= 0 && (edges[inter.Opposite[i]].State == EdgeState.Filled || inter.Opposite[i] == extraEdge))
                    {
                        nextE1 = inter.Edges[i];
                        nextE2 = inter.Opposite[i];
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
            return false;
        }

        private bool IsStraightAssuming(int interIndex, int extraEdge)
        {
            int a, b;
            return IsStraightAssuming(interIndex, extraEdge, out a, out b);
        }

        public void FullClear()
        {
            Clear();
            foreach (Intersection inter in intersections)
            {
                inter.Type = IntersType.Unknown;
            }
            successLookup.Clear();

        }
        public void Clear()
        {
            satisifiedIntersCount = intersections.Count;
            foreach (Edge edge in edges)
            {
                edge.State = EdgeState.Empty;
                edge.Color = 0;
                edge.EdgeSet = 0;
            }
            edgeSets.Clear();
            colorSets.Clear();
            foreach (Cell cell in cells)
            {
                cell.Color = 0;
                cell.FilledCount = 0;
                cell.ExcludedCount = 0;
            }
            cellColorSets.Clear();
            edgeChains.Clear();
            foreach (Intersection inters in intersections)
            {
                inters.ExcludedCount = 0;
                inters.FilledCount = 0;
            }
            edgePairRestrictions = new EdgePairRestriction[edges.Count, edges.Count];
        }

        private int satisifiedIntersCount;
        private List<List<int>> edgeSets = new List<List<int>>();
        private List<List<int>> colorSets = new List<List<int>>();
        private List<List<int>> cellColorSets = new List<List<int>>();
        private List<Chain> edgeChains = new List<Chain>();
        public EdgePairRestriction[] GetEdgePairRestrictionsForEdge(int edge)
        {
            EdgePairRestriction[] restricts = new EdgePairRestriction[edges.Count];
            if (edgePairRestrictions != null)
            {
                for (int i = 0; i < edges.Count; i++)
                {
                    restricts[i] = edgePairRestrictions[edge, i];
                }
            }
            return restricts;
        }
        private EdgePairRestriction[,]? edgePairRestrictions;


        public event EventHandler? PrunedCountProgress;

        public bool AbortPrune = false;

        private bool pruning = false;

        private int[] CreateShuffledTrials(Random rnd)
        {
            int[] trials = new int[intersections.Count];
            for (int i = 0; i < trials.Length; i++)
                trials[i] = i;
            rnd.Shuffle(trials);
            return trials;
        }

        private void TryPruneIntersection(int trial)
        {
            Intersection inter = intersections[trial];
            IntersType oldVal = inter.Type;
            inter.Type = IntersType.Unknown;
            if (TrySolve() != SolveState.Solved)
            {
                inter.Type = oldVal;
            }
            else
            {
                finalSolution = solutionsFound != null && solutionsFound.Count > 0 ? solutionsFound[0] : null;
                finalDepthPatern = solutionDepthPatern != null && solutionDepthPatern.Count > 0 ? solutionDepthPatern[0] : null;
            }
        }

        private bool TryPruneCounts()
        {
            SolveState state = TrySolve();
            if (state != SolveState.Solved)
                return false;
            try
            {
                pruning = true;
                finalSolution = solutionsFound != null && solutionsFound.Count > 0 ? solutionsFound[0] : null;
                finalDepthPatern = solutionDepthPatern != null && solutionDepthPatern.Count > 0 ? solutionDepthPatern[0] : null;
                int[] trials = CreateShuffledTrials(Random.Shared);
                foreach (int trial in trials)
                {
                    if (AbortPrune)
                        break;
                    TryPruneIntersection(trial);
                    PrunedCountProgress?.Invoke(this, EventArgs.Empty);
                }
                return true;
            }
            finally
            {
                pruning = false;
            }
        }


        public void SetRatingCodeOptions(string code)
        {
            IterativeSolverDepth = int.MaxValue;
            IterativeRecMaxDepth = 1;
            ConsiderMultipleLoops = true;
            SolverMethod = SolverMethod.Iterative;
            UseColoring = false;
            UseCellPairs = false;
            UseCellPairsTopLevel = false;
            UseEdgeRestricts = false;
            UseDerivedColoring = false;
            UseCellColoring = false;
            ColoringCheats = false;
            UseMerging = false;
            int state = -1;
            string numAccumulator = string.Empty;
            foreach (char c in code)
            {
                if (char.IsDigit(c) || (numAccumulator.Length == 0 && c == '-'))
                {
                    numAccumulator += c;
                }
                else
                {
                    numAccumulator = HandlePotentialNumber(numAccumulator, state);
                }
                switch (c)
                {
                    case 'F':
                        SolverMethod = SolverMethod.Recursive;
                        state = 1;
                        break;
                    case 'S':
                        SolverMethod = SolverMethod.Iterative;
                        state = 2;
                        break;
                    case 'R':
                        IterativeRecMaxDepth = 2;
                        state = 3;
                        break;
                    case 'C':
                        UseColoring = true;
                        state = 5;
                        break;
                    case 'O':
                        UseCellColoring = true;
                        state = 6;
                        break;
                    case 'N':
                        ConsiderMultipleLoops = false;
                        state = 7;
                        break;
                    case 'H':
                        ColoringCheats = true;
                        state = 8;
                        break;
                    case 'M':
                        UseMerging = true;
                        state = 9;
                        break;
                    case '+':
                        if (state == 5)
                            UseDerivedColoring = true;
                        else if (state == 11)
                            UseCellPairs = true;
                        break;
                    case 'E':
                        UseEdgeRestricts = true;
                        state = 10;
                        break;
                    case 'P':
                        UseCellPairsTopLevel = true;
                        state = 11;
                        break;
                }
            }
            HandlePotentialNumber(numAccumulator, state);
        }

        private string HandlePotentialNumber(string numAccumulator, int state)
        {
            if (string.IsNullOrEmpty(numAccumulator))
                return numAccumulator;
            int value = int.Parse(numAccumulator);
            if (state == 2)
                iterativeSolverDepth = value;
            return string.Empty;
        }

        public bool TypeSatisfied(int i)
        {
            return intersections[i].Type == IntersType.Unknown || intersections[i].Type == DetermineType(i);
        }

        public async Task GenerateAsync(IProgress<int>? progress = null, CancellationToken cancellationToken = default)
        {
            AbortPrune = false;
            bool done = false;
            List<IAction> backup = [];
            Random rnd = Random.Shared;
            while (!done && !cancellationToken.IsCancellationRequested && !AbortPrune)
            {
                try
                {
                    CreatingLoop = true;
                    await GenerateInitialLoopAsync(rnd, backup, cancellationToken);
                    if (cancellationToken.IsCancellationRequested || AbortPrune) break;

                    for (int i = 0; i < 8; i++)
                    {
                        if (cancellationToken.IsCancellationRequested || AbortPrune) break;
                        if (!StepLoopOptimization(i))
                            break;
                    }
                }
                finally
                {
                    CreatingLoop = false;
                }
                Clear();
                if (cancellationToken.IsCancellationRequested || AbortPrune) break;

                done = await TryPruneCountsAsync(progress, cancellationToken);
                if (!done && !cancellationToken.IsCancellationRequested && !AbortPrune)
                {
                    progress?.Report(0);
                    await Task.Delay(1, cancellationToken);
                }
            }
        }

        private async Task<bool> TryPruneCountsAsync(IProgress<int>? progress, CancellationToken cancellationToken)
        {
            if (AbortPrune || cancellationToken.IsCancellationRequested)
                return false;

            await Task.Delay(1, cancellationToken);

            SolveState state = TrySolve();
            if (AbortPrune || cancellationToken.IsCancellationRequested)
                return false;

            if (state != SolveState.Solved)
                return false;

            try
            {
                pruning = true;
                finalSolution = solutionsFound != null && solutionsFound.Count > 0 ? solutionsFound[0] : null;
                finalDepthPatern = solutionDepthPatern != null && solutionDepthPatern.Count > 0 ? solutionDepthPatern[0] : null;
                int[] trials = CreateShuffledTrials(Random.Shared);
                int prunedCount = 0;
                foreach (int trial in trials)
                {
                    if (AbortPrune || cancellationToken.IsCancellationRequested)
                        break;
                    TryPruneIntersection(trial);
                    prunedCount++;
                    progress?.Report(prunedCount);
                    if (prunedCount % 5 == 0)
                    {
                        await Task.Delay(1, cancellationToken);
                    }
                }
                return true;
            }
            finally
            {
                pruning = false;
            }
        }
    }
}