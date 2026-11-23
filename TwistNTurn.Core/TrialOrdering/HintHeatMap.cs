using System;
using System.Collections.Generic;

namespace TwistNTurn.TrialOrdering
{
    /// <summary>
    /// Computes per-edge normalized [0, 1] heat scores for an interactive Masyu puzzle board
    /// using the embedded IncidenceGNN model combined with immediate Masyu local deduction detection.
    /// Higher heat indicates edges that are most promising to examine next.
    /// Not thread-safe: hold one instance per UI session.
    /// </summary>
    public sealed class HintHeatMapScorer
    {
        private readonly IncidenceGnnPolicy _policy;
        private float[] _heat = Array.Empty<float>();
        private bool[] _immediate = Array.Empty<bool>();

        public HintHeatMapScorer(IncidenceGnnModel? model = null)
        {
            _policy = new IncidenceGnnPolicy(model ?? IncidenceGnnModel.Default);
        }

        /// <summary>
        /// Scores every edge in <paramref name="mesh"/> and returns an array of length
        /// <c>mesh.Edges.Count</c> with normalized heat values in [0, 1].
        /// Non-empty edges always receive 0.
        /// </summary>
        public float[] ComputeEdgeHeat(Mesh mesh)
        {
            ArgumentNullException.ThrowIfNull(mesh);

            List<Edge> edges = mesh.Edges;
            int eCount = edges.Count;
            if (_heat.Length != eCount)
            {
                _heat = new float[eCount];
                _immediate = new bool[eCount];
            }
            else
            {
                Array.Clear(_heat, 0, eCount);
                Array.Clear(_immediate, 0, eCount);
            }

            if (eCount == 0)
                return _heat;

            int emptyCount = 0;
            for (int e = 0; e < eCount; e++)
            {
                if (edges[e].State == EdgeState.Empty)
                    emptyCount++;
            }

            if (emptyCount == 0)
                return _heat;

            bool anyImmediate = MarkImmediateLocalDeductions(mesh, _immediate);

            // Temporarily set global solver-mode features to match visible human board state
            // (no internal solver coloring/restriction tables populated) without mutating puzzle state.
            int prevDepth = mesh.IterativeSolverDepth;
            bool prevMulti = mesh.ConsiderMultipleLoops;
            bool prevColor = mesh.UseColoring;
            bool prevCellColor = mesh.UseCellColoring;
            bool prevRestricts = mesh.UseEdgeRestricts;
            bool prevMerging = mesh.UseMerging;
            bool prevCellPairs = mesh.UseCellPairs;
            bool prevCellPairsTop = mesh.UseCellPairsTopLevel;

            try
            {
                mesh.IterativeSolverDepth = 1;
                mesh.ConsiderMultipleLoops = true;
                mesh.UseColoring = false;
                mesh.UseCellColoring = false;
                mesh.UseEdgeRestricts = false;
                mesh.UseMerging = false;
                mesh.UseCellPairs = false;
                mesh.UseCellPairsTopLevel = false;

                _policy.ScoreAllNodes(mesh, iterativeRecDepth: 1);
            }
            finally
            {
                mesh.IterativeSolverDepth = prevDepth;
                mesh.ConsiderMultipleLoops = prevMulti;
                mesh.UseColoring = prevColor;
                mesh.UseCellColoring = prevCellColor;
                mesh.UseEdgeRestricts = prevRestricts;
                mesh.UseMerging = prevMerging;
                mesh.UseCellPairs = prevCellPairs;
                mesh.UseCellPairsTopLevel = prevCellPairsTop;
            }

            float[] rawScores = _policy.EdgeScores;

            // Compute mean, std, and max over all empty edges
            double sum = 0.0;
            float maxRaw = float.NegativeInfinity;
            for (int e = 0; e < eCount; e++)
            {
                if (edges[e].State != EdgeState.Empty)
                    continue;
                float s = rawScores[e];
                if (!float.IsFinite(s))
                    continue;
                sum += s;
                if (s > maxRaw)
                    maxRaw = s;
            }

            float mean = (float)(sum / emptyCount);
            double varSum = 0.0;
            for (int e = 0; e < eCount; e++)
            {
                if (edges[e].State != EdgeState.Empty)
                    continue;
                float s = rawScores[e];
                if (!float.IsFinite(s))
                    continue;
                double d = s - mean;
                varSum += d * d;
            }

            float std = MathF.Sqrt((float)(varSum / emptyCount));
            float zMax = std > 1e-4f ? (maxRaw - mean) / std : 0f;
            float spreadScale = zMax > 0f ? MathF.Min(1.0f, zMax / 1.5f) : 0f;
            float zDenom = MathF.Max(zMax, 1.5f);

            // When immediate 0-step local deductions exist, reserve the top heat band [0.90, 1.0]
            // for them while keeping strategic GNN hotspots visible in [0, 0.72].
            float gnnCap = anyImmediate ? 0.72f : 1.0f;

            for (int e = 0; e < eCount; e++)
            {
                if (edges[e].State != EdgeState.Empty)
                    continue;

                float gnnHeat = 0f;
                if (std > 1e-4f && spreadScale > 0f)
                {
                    float s = rawScores[e];
                    if (float.IsFinite(s) && s > mean)
                    {
                        float z = (s - mean) / std;
                        float t = Math.Clamp(z / zDenom, 0f, 1f);
                        gnnHeat = spreadScale * MathF.Pow(t, 1.25f);
                    }
                }

                if (_immediate[e])
                {
                    _heat[e] = 0.90f + 0.10f * gnnHeat;
                }
                else
                {
                    _heat[e] = gnnHeat * gnnCap;
                }
            }

            return _heat;
        }

        private static bool MarkImmediateLocalDeductions(Mesh mesh, bool[] immediate)
        {
            List<Edge> edges = mesh.Edges;
            List<Intersection> inters = mesh.Intersections;
            bool any = false;

            for (int v = 0; v < inters.Count; v++)
            {
                Intersection inter = inters[v];
                List<int> vEdges = inter.Edges;
                List<int> vOpp = inter.Opposite;
                int deg = vEdges.Count;
                int rem = deg - inter.FilledCount - inter.ExcludedCount;

                // 1. Standard loop degree rules
                bool degreeForced =
                    (inter.FilledCount == 2 && rem > 0) ||
                    (inter.FilledCount == 0 && rem == 1) ||
                    ((inter.FilledCount > 0 || inter.Type != IntersType.Unknown) && (deg - inter.ExcludedCount == 2) && inter.FilledCount < 2);

                if (degreeForced)
                {
                    for (int i = 0; i < deg; i++)
                    {
                        int eIdx = vEdges[i];
                        if (edges[eIdx].State == EdgeState.Empty)
                        {
                            immediate[eIdx] = true;
                            any = true;
                        }
                    }
                }

                // 2. Black pearl (IntersType.Turn): must turn at pearl and extend straight for 2 edges on both legs
                if (inter.Type == IntersType.Turn)
                {
                    for (int i = 0; i < deg; i++)
                    {
                        int eIdx = vEdges[i];
                        Edge edge = edges[eIdx];
                        int opp = i < vOpp.Count ? vOpp[i] : -1;

                        if (edge.State == EdgeState.Filled)
                        {
                            // Opposite edge at black pearl must be excluded
                            if (opp >= 0 && edges[opp].State == EdgeState.Empty)
                            {
                                immediate[opp] = true;
                                any = true;
                            }

                            // Straight-through continuation at the far end must be filled
                            int otherV = GetOtherVertex(edge, v);
                            if (otherV >= 0)
                            {
                                int cont = GetOppositeAtVertex(inters[otherV], eIdx);
                                if (cont >= 0 && edges[cont].State == EdgeState.Empty)
                                {
                                    immediate[cont] = true;
                                    any = true;
                                }
                            }
                        }
                        else if (edge.State == EdgeState.Empty)
                        {
                            // Opposite is missing or excluded on a degree-3/4 vertex -> must fill this leg
                            if (((inter.ExcludedCount == 1 && deg == 4) || deg == 3)
                                && (opp == -1 || edges[opp].State == EdgeState.Excluded))
                            {
                                immediate[eIdx] = true;
                                any = true;
                                continue;
                            }

                            int otherV = GetOtherVertex(edge, v);
                            if (otherV >= 0)
                            {
                                Intersection otherInter = inters[otherV];
                                int cont = GetOppositeAtVertex(otherInter, eIdx);

                                // Far end cannot continue straight -> this edge must be excluded
                                if (cont == -1 || edges[cont].State == EdgeState.Excluded)
                                {
                                    immediate[eIdx] = true;
                                    any = true;
                                    continue;
                                }

                                // Far end already has a perpendicular filled edge -> turning too early, must exclude
                                if (otherInter.FilledCount == 1)
                                {
                                    for (int j = 0; j < otherInter.Edges.Count; j++)
                                    {
                                        int oe = otherInter.Edges[j];
                                        if (oe != eIdx && oe != cont && edges[oe].State == EdgeState.Filled)
                                        {
                                            immediate[eIdx] = true;
                                            any = true;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                // 3. White pearl (IntersType.Twist): must go straight through pearl and turn on at least one side
                else if (inter.Type == IntersType.Twist && deg > 2)
                {
                    for (int i = 0; i < deg; i++)
                    {
                        int eIdx = vEdges[i];
                        Edge edge = edges[eIdx];
                        int opp = i < vOpp.Count ? vOpp[i] : -1;

                        if (edge.State == EdgeState.Filled)
                        {
                            // Opposite edge at white pearl must be filled
                            if (opp >= 0 && edges[opp].State == EdgeState.Empty)
                            {
                                immediate[opp] = true;
                                any = true;
                            }

                            // If the path continues straight past one side, the opposite side must turn
                            int otherV = GetOtherVertex(edge, v);
                            if (otherV >= 0 && opp >= 0)
                            {
                                int cont = GetOppositeAtVertex(inters[otherV], eIdx);
                                if (cont >= 0 && edges[cont].State == EdgeState.Filled)
                                {
                                    int finalV = GetOtherVertex(edges[opp], v);
                                    if (finalV >= 0)
                                    {
                                        int finalCont = GetOppositeAtVertex(inters[finalV], opp);
                                        if (finalCont >= 0 && edges[finalCont].State == EdgeState.Empty)
                                        {
                                            immediate[finalCont] = true;
                                            any = true;
                                        }
                                    }
                                }
                            }
                        }
                        else if (edge.State == EdgeState.Empty)
                        {
                            // Missing or excluded opposite -> cannot go straight along this axis, must exclude
                            if (opp == -1 || edges[opp].State == EdgeState.Excluded)
                            {
                                immediate[eIdx] = true;
                                any = true;
                            }
                        }
                    }
                }
            }

            return any;
        }

        private static int GetOtherVertex(Edge edge, int vIdx)
        {
            if (edge.Intersections == null || edge.Intersections.Length < 2)
                return -1;
            return edge.Intersections[0] == vIdx ? edge.Intersections[1] : edge.Intersections[0];
        }

        private static int GetOppositeAtVertex(Intersection inter, int edgeIdx)
        {
            List<int> vEdges = inter.Edges;
            List<int> vOpp = inter.Opposite;
            for (int i = 0; i < vEdges.Count; i++)
            {
                if (vEdges[i] == edgeIdx)
                    return i < vOpp.Count ? vOpp[i] : -1;
            }
            return -1;
        }
    }
}
