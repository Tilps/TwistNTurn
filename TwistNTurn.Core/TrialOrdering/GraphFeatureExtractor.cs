using System;
using System.Collections.Generic;

namespace TwistNTurn.TrialOrdering
{
    /// <summary>
    /// Cached static planar incidence and collinear-continuation topology for a <see cref="Mesh"/>,
    /// avoiding repeated list traversals and allocations during trial scoring.
    /// Precomputes straight-through opposite edges (<c>EdgeOpp0</c>/<c>EdgeOpp1</c>) in O(E).
    /// </summary>
    public sealed class MeshTopologyCache
    {
        public readonly int NumVerts;
        public readonly int NumEdges;
        public readonly int NumCells;

        public readonly int[] EdgeV0;
        public readonly int[] EdgeV1;
        public readonly int[] EdgeC0;
        public readonly int[] EdgeC1;

        // Collinear straight-through neighbor edge at v0 and v1 (-1 if none), and its far intersection
        public readonly int[] EdgeOpp0;
        public readonly int[] EdgeOpp1;
        public readonly int[] EdgeOppV0;
        public readonly int[] EdgeOppV1;

        public readonly int[] VertEdgeOff;
        public readonly int[] VertEdgeIdx;
        public readonly int[] VertEdgeOpp;

        public readonly int[] CellEdgeOff;
        public readonly int[] CellEdgeIdx;
        public readonly int[] CellVertOff;
        public readonly int[] CellVertIdx;

        public readonly int[] EdgeNbrOff;
        public readonly int[] EdgeNbrIdx;

        public MeshTopologyCache(Mesh mesh)
        {
            NumVerts = mesh.Intersections.Count;
            NumEdges = mesh.Edges.Count;
            NumCells = mesh.Cells.Count;

            EdgeV0 = new int[NumEdges];
            EdgeV1 = new int[NumEdges];
            EdgeC0 = new int[NumEdges];
            EdgeC1 = new int[NumEdges];
            EdgeOpp0 = new int[NumEdges];
            EdgeOpp1 = new int[NumEdges];
            EdgeOppV0 = new int[NumEdges];
            EdgeOppV1 = new int[NumEdges];

            for (int e = 0; e < NumEdges; e++)
            {
                Edge edge = mesh.Edges[e];
                EdgeV0[e] = edge.Intersections[0];
                EdgeV1[e] = edge.Intersections[1];
                EdgeC0[e] = edge.Cells.Count > 0 ? edge.Cells[0] : -1;
                EdgeC1[e] = edge.Cells.Count > 1 ? edge.Cells[1] : -1;
                EdgeOpp0[e] = -1;
                EdgeOpp1[e] = -1;
                EdgeOppV0[e] = -1;
                EdgeOppV1[e] = -1;
            }

            VertEdgeOff = new int[NumVerts + 1];
            int totalVE = 0;
            for (int v = 0; v < NumVerts; v++)
            {
                VertEdgeOff[v] = totalVE;
                totalVE += mesh.Intersections[v].Edges.Count;
            }
            VertEdgeOff[NumVerts] = totalVE;
            VertEdgeIdx = new int[totalVE];
            VertEdgeOpp = new int[totalVE];
            for (int v = 0; v < NumVerts; v++)
            {
                int off = VertEdgeOff[v];
                Intersection inter = mesh.Intersections[v];
                List<int> ve = inter.Edges;
                List<int> vopp = inter.Opposite;
                for (int i = 0; i < ve.Count; i++)
                {
                    int eIdx = ve[i];
                    int opp = i < vopp.Count ? vopp[i] : -1;
                    VertEdgeIdx[off + i] = eIdx;
                    VertEdgeOpp[off + i] = opp;
                    if (EdgeV0[eIdx] == v)
                    {
                        EdgeOpp0[eIdx] = opp;
                        if (opp >= 0)
                            EdgeOppV0[eIdx] = EdgeV0[opp] == v ? EdgeV1[opp] : EdgeV0[opp];
                    }
                    else
                    {
                        EdgeOpp1[eIdx] = opp;
                        if (opp >= 0)
                            EdgeOppV1[eIdx] = EdgeV0[opp] == v ? EdgeV1[opp] : EdgeV0[opp];
                    }
                }
            }

            CellEdgeOff = new int[NumCells + 1];
            CellVertOff = new int[NumCells + 1];
            int totalCE = 0, totalCV = 0;
            for (int c = 0; c < NumCells; c++)
            {
                CellEdgeOff[c] = totalCE;
                totalCE += mesh.Cells[c].Edges.Count;
                CellVertOff[c] = totalCV;
                totalCV += mesh.Cells[c].Intersections.Count;
            }
            CellEdgeOff[NumCells] = totalCE;
            CellVertOff[NumCells] = totalCV;
            CellEdgeIdx = new int[totalCE];
            CellVertIdx = new int[totalCV];
            for (int c = 0; c < NumCells; c++)
            {
                int offE = CellEdgeOff[c];
                List<int> ce = mesh.Cells[c].Edges;
                for (int i = 0; i < ce.Count; i++) CellEdgeIdx[offE + i] = ce[i];
                int offV = CellVertOff[c];
                List<int> cv = mesh.Cells[c].Intersections;
                for (int i = 0; i < cv.Count; i++) CellVertIdx[offV + i] = cv[i];
            }

            List<int> nbrList = new List<int>(NumEdges * 8);
            EdgeNbrOff = new int[NumEdges + 1];
            int[] seenMark = new int[NumEdges];
            for (int e = 0; e < NumEdges; e++)
            {
                EdgeNbrOff[e] = nbrList.Count;
                int mark = e + 1;
                seenMark[e] = mark;

                int v0 = EdgeV0[e], v1 = EdgeV1[e];
                for (int k = VertEdgeOff[v0]; k < VertEdgeOff[v0 + 1]; k++)
                {
                    int oe = VertEdgeIdx[k];
                    if (seenMark[oe] != mark) { seenMark[oe] = mark; nbrList.Add(oe); }
                }
                for (int k = VertEdgeOff[v1]; k < VertEdgeOff[v1 + 1]; k++)
                {
                    int oe = VertEdgeIdx[k];
                    if (seenMark[oe] != mark) { seenMark[oe] = mark; nbrList.Add(oe); }
                }
                int c0 = EdgeC0[e], c1 = EdgeC1[e];
                if (c0 >= 0)
                {
                    for (int k = CellEdgeOff[c0]; k < CellEdgeOff[c0 + 1]; k++)
                    {
                        int oe = CellEdgeIdx[k];
                        if (seenMark[oe] != mark) { seenMark[oe] = mark; nbrList.Add(oe); }
                    }
                }
                if (c1 >= 0)
                {
                    for (int k = CellEdgeOff[c1]; k < CellEdgeOff[c1 + 1]; k++)
                    {
                        int oe = CellEdgeIdx[k];
                        if (seenMark[oe] != mark) { seenMark[oe] = mark; nbrList.Add(oe); }
                    }
                }
            }
            EdgeNbrOff[NumEdges] = nbrList.Count;
            EdgeNbrIdx = nbrList.ToArray();
        }
    }

    /// <summary>
    /// Zero-allocation O(V + E + C) feature extractor capturing Masyu white/black pearl clues,
    /// collinear straight-continuation geometry, solver-strength flags, internal coloring sets,
    /// loop edge sets, and local O(E) projections of <c>edgePairRestrictions</c>.
    /// </summary>
    public static class GraphFeatureExtractor
    {
        public const int GlobalFeatDim = 12;
        public const int EdgeFeatDim = 24;
        public const int CellFeatDim = 9;
        public const int VertFeatDim = 14;

        private static readonly float[] Log2Div5 = BuildLogTable(512, 5.0f);
        private static readonly float[] Log2Div6 = BuildLogTable(512, 6.0f);

        private static float[] BuildLogTable(int n, float div)
        {
            float[] t = new float[n];
            for (int i = 0; i < n; i++) t[i] = MathF.Log2(i + 1) / div;
            return t;
        }

        private static float FastLog2Div5(int n) => (uint)n < (uint)Log2Div5.Length ? Log2Div5[n] : MathF.Log2(n + 1) / 5.0f;
        private static float FastLog2Div6(int n) => (uint)n < (uint)Log2Div6.Length ? Log2Div6[n] : MathF.Log2(n + 1) / 6.0f;

        public static void Extract(
            Mesh mesh,
            MeshTopologyCache topo,
            float[] globalFeats,
            float[] edgeFeats,
            float[] cellFeats,
            float[] vertFeats,
            int[] edgeColors,
            int[] cellColors,
            int[] vertEdgeSetScratch,
            int iterativeRecDepth)
        {
            int V = topo.NumVerts;
            int E = topo.NumEdges;
            int C = topo.NumCells;

            List<Edge> edges = mesh.Edges;
            List<Cell> cells = mesh.Cells;
            List<Intersection> inters = mesh.Intersections;
            IReadOnlyList<List<int>> edgeSets = mesh.EdgeSets;
            IReadOnlyList<List<int>> colorSets = mesh.ColorSets;
            IReadOnlyList<List<int>> cellColorSets = mesh.CellColorSets;

            int emptyEdges = 0;
            for (int e = 0; e < E; e++)
                if (edges[e].State == EdgeState.Empty) emptyEdges++;

            int activeEdgeSets = 0;
            for (int i = 0; i < edgeSets.Count; i++)
                if (edgeSets[i].Count > 0) activeEdgeSets++;

            int totalPearls = 0;
            int unsatisfiedPearls = 0;

            for (int v = 0; v < V; v++)
            {
                Intersection inter = inters[v];
                int deg = inter.Edges.Count;
                int rem = deg - inter.FilledCount - inter.ExcludedCount;
                int es = 0;
                int start = topo.VertEdgeOff[v], end = topo.VertEdgeOff[v + 1];
                int straightPairsTotal = 0;
                int straightPairsOpen = 0;
                for (int k = start; k < end; k++)
                {
                    int eIdx = topo.VertEdgeIdx[k];
                    Edge ie = edges[eIdx];
                    if (es == 0 && ie.State == EdgeState.Filled)
                        es = ie.EdgeSet;
                    int opp = topo.VertEdgeOpp[k];
                    if (opp > eIdx)
                    {
                        straightPairsTotal++;
                        if (ie.State != EdgeState.Excluded && edges[opp].State != EdgeState.Excluded)
                            straightPairsOpen++;
                    }
                }
                vertEdgeSetScratch[v] = es;

                bool isTwist = inter.Type == IntersType.Twist;
                bool isTurn = inter.Type == IntersType.Turn;
                bool hasPearl = isTwist || isTurn;
                bool satisfied = false;
                if (hasPearl)
                {
                    totalPearls++;
                    satisfied = inter.FilledCount == 2 && mesh.TypeSatisfied(v);
                    if (!satisfied) unsatisfiedPearls++;
                }

                int o = v * VertFeatDim;
                vertFeats[o + 0] = deg / 6.0f;
                vertFeats[o + 1] = inter.FilledCount == 0 ? 1f : 0f;
                vertFeats[o + 2] = inter.FilledCount == 1 ? 1f : 0f;
                vertFeats[o + 3] = inter.FilledCount == 2 ? 1f : 0f;
                vertFeats[o + 4] = deg > 0 ? (float)rem / deg : 0f;
                vertFeats[o + 5] = rem == 1 ? 1f : 0f;
                vertFeats[o + 6] = rem == 2 ? 1f : 0f;
                vertFeats[o + 7] = es != 0 ? 1f : 0f;
                vertFeats[o + 8] = isTwist ? 1f : 0f;
                vertFeats[o + 9] = isTurn ? 1f : 0f;
                vertFeats[o + 10] = (hasPearl && satisfied) ? 1f : 0f;
                vertFeats[o + 11] = (hasPearl && inter.FilledCount == 0) ? 1f : 0f;
                vertFeats[o + 12] = (hasPearl && inter.FilledCount == 1) ? 1f : 0f;
                vertFeats[o + 13] = straightPairsTotal > 0 ? (float)straightPairsOpen / straightPairsTotal : 0f;
            }

            globalFeats[0] = mesh.IterativeSolverDepth > 1000 ? 1.0f : MathF.Log2(Math.Max(0, mesh.IterativeSolverDepth) + 1) / 8.0f;
            globalFeats[1] = mesh.IterativeSolverDepth >= E ? 1.0f : 0.0f;
            globalFeats[2] = iterativeRecDepth >= 2 ? 1.0f : 0.0f;
            globalFeats[3] = mesh.ConsiderMultipleLoops ? 1.0f : 0.0f;
            globalFeats[4] = mesh.UseColoring ? 1.0f : 0.0f;
            globalFeats[5] = mesh.UseCellColoring ? 1.0f : 0.0f;
            globalFeats[6] = mesh.UseEdgeRestricts ? 1.0f : 0.0f;
            globalFeats[7] = mesh.UseMerging ? 1.0f : 0.0f;
            globalFeats[8] = (mesh.UseCellPairs || mesh.UseCellPairsTopLevel) ? 1.0f : 0.0f;
            globalFeats[9] = E > 0 ? (float)emptyEdges / E : 0f;
            globalFeats[10] = totalPearls > 0 ? (float)unsatisfiedPearls / totalPearls : 0f;
            globalFeats[11] = MathF.Min(activeEdgeSets, 10) / 10.0f;

            bool useCellColoring = mesh.UseCellColoring;
            for (int c = 0; c < C; c++)
            {
                Cell cell = cells[c];
                int deg = cell.Edges.Count;
                int rem = deg - cell.FilledCount - cell.ExcludedCount;
                int col = useCellColoring ? cell.Color : 0;
                cellColors[c] = col;

                int adjColoredEmpty = 0;
                int emptyBoundary = 0;
                int cStart = topo.CellEdgeOff[c], cEnd = topo.CellEdgeOff[c + 1];
                for (int k = cStart; k < cEnd; k++)
                {
                    int eIdx = topo.CellEdgeIdx[k];
                    if (edges[eIdx].State == EdgeState.Empty)
                    {
                        emptyBoundary++;
                        int oc = topo.EdgeC0[eIdx] == c ? topo.EdgeC1[eIdx] : topo.EdgeC0[eIdx];
                        if (oc == -1 || (useCellColoring && cells[oc].Color != 0))
                            adjColoredEmpty++;
                    }
                }

                int pearlCorners = 0;
                int cvStart = topo.CellVertOff[c], cvEnd = topo.CellVertOff[c + 1];
                int cvCnt = cvEnd - cvStart;
                for (int k = cvStart; k < cvEnd; k++)
                {
                    if (inters[topo.CellVertIdx[k]].Type != IntersType.Unknown)
                        pearlCorners++;
                }

                int o = c * CellFeatDim;
                cellFeats[o + 0] = deg / 8.0f;
                cellFeats[o + 1] = deg > 0 ? (float)cell.FilledCount / deg : 0f;
                cellFeats[o + 2] = deg > 0 ? (float)cell.ExcludedCount / deg : 0f;
                cellFeats[o + 3] = deg > 0 ? (float)rem / deg : 0f;
                cellFeats[o + 4] = Math.Abs(col) == 1 ? 1f : 0f;
                cellFeats[o + 5] = Math.Abs(col) > 1 ? 1f : 0f;
                int cellColorIdx = Math.Abs(col) - 1;
                cellFeats[o + 6] = (cellColorIdx >= 0 && cellColorIdx < cellColorSets.Count)
                    ? FastLog2Div5(cellColorSets[cellColorIdx].Count)
                    : 0f;
                cellFeats[o + 7] = emptyBoundary > 0 ? (float)adjColoredEmpty / emptyBoundary : 0f;
                cellFeats[o + 8] = cvCnt > 0 ? (float)pearlCorners / cvCnt : 0f;
            }

            bool useRestricts = mesh.UseEdgeRestricts;
            bool useColoring = mesh.UseColoring;
            for (int e = 0; e < E; e++)
            {
                Edge edge = edges[e];
                int col = useColoring ? edge.Color : 0;
                edgeColors[e] = col;

                int sameCol = 0, oppCol = 0, notBoth = 0, notNeither = 0;
                int nStart = topo.EdgeNbrOff[e], nEnd = topo.EdgeNbrOff[e + 1];
                int nCnt = nEnd - nStart;
                if (useColoring || useRestricts)
                {
                    for (int k = nStart; k < nEnd; k++)
                    {
                        int oe = topo.EdgeNbrIdx[k];
                        if (useColoring && col != 0)
                        {
                            int ocol = edges[oe].Color;
                            if (ocol == col) sameCol++;
                            else if (ocol == -col) oppCol++;
                        }
                        if (useRestricts)
                        {
                            EdgePairRestriction r = mesh.GetEdgePairRestriction(e, oe);
                            if (r == EdgePairRestriction.NotBoth) notBoth++;
                            else if (r == EdgePairRestriction.NotNeither) notNeither++;
                        }
                    }
                }

                int v0 = topo.EdgeV0[e], v1 = topo.EdgeV1[e];
                int es0 = vertEdgeSetScratch[v0];
                int es1 = vertEdgeSetScratch[v1];
                int c0 = topo.EdgeC0[e], c1 = topo.EdgeC1[e];
                int c0Col = c0 >= 0 ? (useCellColoring ? cells[c0].Color : 0) : 1;
                int c1Col = c1 >= 0 ? (useCellColoring ? cells[c1].Color : 0) : 1;

                bool cellsSame = c0Col != 0 && c0Col == c1Col;
                bool cellsOpp = c0Col != 0 && c0Col == -c1Col;
                bool cellsDiffNonZero = c0Col != 0 && c1Col != 0 && Math.Abs(c0Col) != Math.Abs(c1Col);

                int opp0 = topo.EdgeOpp0[e], opp1 = topo.EdgeOpp1[e];
                int oppAvail = (opp0 >= 0 ? 1 : 0) + (opp1 >= 0 ? 1 : 0);
                int oppFilled = (opp0 >= 0 && edges[opp0].State == EdgeState.Filled ? 1 : 0)
                              + (opp1 >= 0 && edges[opp1].State == EdgeState.Filled ? 1 : 0);
                int oppExcl = (opp0 >= 0 && edges[opp0].State == EdgeState.Excluded ? 1 : 0)
                            + (opp1 >= 0 && edges[opp1].State == EdgeState.Excluded ? 1 : 0);
                IntersType t0 = inters[v0].Type, t1 = inters[v1].Type;
                int ov0 = topo.EdgeOppV0[e], ov1 = topo.EdgeOppV1[e];
                bool twoHopBlack = (ov0 >= 0 && inters[ov0].Type == IntersType.Turn)
                                || (ov1 >= 0 && inters[ov1].Type == IntersType.Turn);

                int o = e * EdgeFeatDim;
                edgeFeats[o + 0] = edge.State == EdgeState.Empty ? 1f : 0f;
                edgeFeats[o + 1] = edge.State == EdgeState.Filled ? 1f : 0f;
                edgeFeats[o + 2] = edge.State == EdgeState.Excluded ? 1f : 0f;
                edgeFeats[o + 3] = c1 < 0 ? 1f : 0f;
                edgeFeats[o + 4] = col != 0 ? 1f : 0f;
                int edgeColorIdx = Math.Abs(col) - 1;
                edgeFeats[o + 5] = (edgeColorIdx >= 0 && edgeColorIdx < colorSets.Count)
                    ? FastLog2Div5(colorSets[edgeColorIdx].Count)
                    : 0f;
                edgeFeats[o + 6] = nCnt > 0 ? (float)sameCol / nCnt : 0f;
                edgeFeats[o + 7] = nCnt > 0 ? (float)oppCol / nCnt : 0f;
                edgeFeats[o + 8] = MathF.Min(notBoth, 6) / 6.0f;
                edgeFeats[o + 9] = MathF.Min(notNeither, 6) / 6.0f;
                edgeFeats[o + 10] = edge.EdgeSet != 0 ? 1f : 0f;
                int edgeSetIdx = edge.EdgeSet - 1;
                edgeFeats[o + 11] = (edgeSetIdx >= 0 && edgeSetIdx < edgeSets.Count)
                    ? FastLog2Div6(edgeSets[edgeSetIdx].Count)
                    : 0f;
                edgeFeats[o + 12] = (es0 != 0 && es0 == es1) ? 1f : 0f;
                edgeFeats[o + 13] = (es0 != 0 && es1 != 0 && es0 != es1) ? 1f : 0f;
                edgeFeats[o + 14] = ((es0 != 0) ^ (es1 != 0)) ? 1f : 0f;
                edgeFeats[o + 15] = cellsSame ? 1f : 0f;
                edgeFeats[o + 16] = cellsOpp ? 1f : 0f;
                edgeFeats[o + 17] = cellsDiffNonZero ? 1f : 0f;
                edgeFeats[o + 18] = oppAvail * 0.5f;
                edgeFeats[o + 19] = oppFilled * 0.5f;
                edgeFeats[o + 20] = oppExcl * 0.5f;
                edgeFeats[o + 21] = (t0 == IntersType.Twist || t1 == IntersType.Twist) ? 1f : 0f;
                edgeFeats[o + 22] = (t0 == IntersType.Turn || t1 == IntersType.Turn) ? 1f : 0f;
                edgeFeats[o + 23] = twoHopBlack ? 1f : 0f;
            }
        }
    }
}
