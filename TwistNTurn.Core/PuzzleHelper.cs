using System;
using System.Collections.Generic;
using System.Globalization;

namespace TwistNTurn
{
    public sealed record PuzzleAnalysisResult(
        string SolveStatus,
        string TechniqueTier,
        int BranchDepth,
        int WhitePearls,
        int BlackPearls,
        int TotalIntersections,
        int EdgeCount,
        int? LookaheadDepth = null)
    {
        public int ClueCount => WhitePearls + BlackPearls;
    }

    public static class PuzzleHelper
    {
        public const int MaxGridDimension = 100;

        public static Mesh MakeMesh(int width, int height, MeshType type, int difficulty)
        {
            Mesh mesh = new Mesh(width, height, type);
            mesh.ConsiderMultipleLoops = difficulty > 0;
            mesh.IterativeRecMaxDepth = 1;
            if (difficulty > 1)
            {
                mesh.UseCellColoring = true;
                mesh.UseCellColoringTrials = true;
                mesh.UseColoring = true;
                mesh.UseEdgeRestricts = true;
                mesh.UseDerivedColoring = true;
                mesh.UseMerging = true;
                mesh.UseCellPairsTopLevel = true;
            }
            if (difficulty < 5)
            {
                mesh.IterativeSolverDepth = Math.Max(0, difficulty - 2);
            }
            else
            {
                mesh.SolverMethod = SolverMethod.Recursive;
            }
            mesh.GenerateBoringFraction = 0.01;
            return mesh;
        }

        public static Mesh CreateFullSolverMesh(Mesh source)
        {
            var solver = new Mesh(source);
            solver.ConsiderMultipleLoops = true;
            solver.IterativeRecMaxDepth = 1;
            solver.SolverMethod = SolverMethod.Recursive;
            solver.UseCellColoring = true;
            solver.UseCellColoringTrials = true;
            solver.UseColoring = true;
            solver.UseEdgeRestricts = true;
            solver.UseDerivedColoring = true;
            solver.UseMerging = true;
            solver.UseCellPairs = true;
            solver.UseCellPairsTopLevel = true;
            return solver;
        }

        public static PuzzleAnalysisResult AnalyzePuzzle(Mesh source)
        {
            int whitePearls = 0;
            int blackPearls = 0;
            foreach (var inter in source.Intersections)
            {
                if (inter.Type == IntersType.Twist) whitePearls++;
                else if (inter.Type == IntersType.Turn) blackPearls++;
            }
            int totalIntersections = source.Intersections.Count;
            int edgeCount = source.Edges.Count;

            Mesh solver = CreateFullSolverMesh(source);
            solver.Clear();

            SolveState fullResult = solver.TrySolve();
            if (fullResult == SolveState.MultipleSolutions)
            {
                return new PuzzleAnalysisResult("Multiple Solutions", "Non-unique", 0, whitePearls, blackPearls, totalIntersections, edgeCount);
            }
            if (fullResult != SolveState.Solved)
            {
                return new PuzzleAnalysisResult("No Solution", "Unsolvable", 0, whitePearls, blackPearls, totalIntersections, edgeCount);
            }

            // RecursiveTrySolve pushes a base frame before checking trials, so branch count is Length - 1.
            int branchDepth = Math.Max(0, (solver.DepthPatern?.Length ?? 0) - 1);
            string tier = $"Recursive Search (F, Depth {branchDepth})";
            int? lookaheadDepth = null;

            solver.Clear();
            solver.SetRatingCodeOptions("S");
            if (solver.TrySolve() == SolveState.Solved)
            {
                int min = -1;
                int max = solver.Edges.Count;
                while (min < max - 1)
                {
                    int mid = (min + max) / 2;
                    solver.Clear();
                    solver.IterativeSolverDepth = mid;
                    if (solver.TrySolve() == SolveState.Solved)
                        max = mid;
                    else
                        min = mid;
                }
                lookaheadDepth = max;
                tier = max == 0 ? "Direct Deductions (S0)" : $"Iterative Lookahead (S{max})";
            }
            else
            {
                (string code, string label)[] tiers =
                [
                    ("SC", "Edge Coloring (SC)"),
                    ("SC+", "Derived Coloring (SC+)"),
                    ("SCE", "Edge Restrictions (SCE)"),
                    ("SCOM", "Cell Coloring & Merging (SCOM)")
                ];
                foreach (var (code, label) in tiers)
                {
                    solver.Clear();
                    solver.SetRatingCodeOptions(code);
                    if (solver.TrySolve() == SolveState.Solved)
                    {
                        tier = label;
                        break;
                    }
                }
            }

            return new PuzzleAnalysisResult("Unique Solution", tier, branchDepth, whitePearls, blackPearls, totalIntersections, edgeCount, lookaheadDepth);
        }

        public static void GetDefaultSize(MeshType type, out int width, out int height)
        {
            (width, height) = type switch
            {
                MeshType.Hexagonal2 => (5, 5),
                MeshType.Hexagonal3 or MeshType.Hexagonal4 => (3, 3),
                _ => (10, 10)
            };
        }

        public static bool ParseSize(string val, MeshType type, out int width, out int height)
        {
            GetDefaultSize(type, out width, out height);
            if (string.IsNullOrWhiteSpace(val))
                return true;

            ReadOnlySpan<char> span = val.AsSpan().Trim();
            int sepIndex = span.IndexOfAny('x', 'X', '*');
            if (sepIndex >= 0)
            {
                ReadOnlySpan<char> wSpan = span[..sepIndex].Trim();
                ReadOnlySpan<char> hSpan = span[(sepIndex + 1)..].Trim();
                if (!int.TryParse(wSpan, NumberStyles.Integer, CultureInfo.InvariantCulture, out width)
                    || !int.TryParse(hSpan, NumberStyles.Integer, CultureInfo.InvariantCulture, out height))
                {
                    return false;
                }
                return width > 0 && height > 0 && width <= MaxGridDimension && height <= MaxGridDimension;
            }
            else
            {
                if (!int.TryParse(span, NumberStyles.Integer, CultureInfo.InvariantCulture, out width))
                    return false;
                height = width;
                return width > 0 && width <= MaxGridDimension;
            }
        }

        public static string StringFromMeshType(MeshType type) => type switch
        {
            MeshType.Square => "Square",
            MeshType.SquareSymmetrical => "Square Symmetrical",
            MeshType.Hexagonal2 => "Hexagon 2",
            MeshType.Hexagonal3 => "Hexagon 3",
            MeshType.Hexagonal4 => "Hexagon 4",
            _ => type.ToString()
        };

        public static MeshType MeshTypeFromString(string typeName)
        {
            if (string.IsNullOrWhiteSpace(typeName))
                return MeshType.Square;
            string normalized = typeName.Trim();
            return normalized.Replace(" ", "").Replace("-", "").Replace("_", "").ToLowerInvariant() switch
            {
                "square" or "sq" => MeshType.Square,
                "squaresymmetrical" or "sqs" => MeshType.SquareSymmetrical,
                "hexagon2" or "hexagonal2" or "hex2" => MeshType.Hexagonal2,
                "hexagon3" or "hexagonal3" or "hex3" => MeshType.Hexagonal3,
                "hexagon4" or "hexagonal4" or "hex4" => MeshType.Hexagonal4,
                _ => Enum.TryParse<MeshType>(normalized, ignoreCase: true, out var parsed) && Enum.IsDefined(parsed)
                    ? parsed
                    : MeshType.Square
            };
        }

        public static bool CheckIsSolved(Mesh mesh)
        {
            if (mesh == null || mesh.Intersections.Count == 0 || mesh.Edges.Count == 0)
                return false;

            mesh.RecalculateCounts();

            for (int i = 0; i < mesh.Intersections.Count; i++)
            {
                if (mesh.Intersections[i].Type != IntersType.Unknown)
                {
                    if (!mesh.TypeSatisfied(i))
                        return false;
                }
            }

            int filledEdgeCount = 0;
            for (int i = 0; i < mesh.Edges.Count; i++)
            {
                if (mesh.Edges[i].State == EdgeState.Filled)
                    filledEdgeCount++;
            }
            if (filledEdgeCount < 3)
                return false;

            for (int i = 0; i < mesh.Intersections.Count; i++)
            {
                int deg = mesh.Intersections[i].FilledCount;
                if (deg != 0 && deg != 2)
                    return false;
            }

            // Verify that all filled edges form a single connected cycle
            int startIntersection = -1;
            for (int i = 0; i < mesh.Intersections.Count; i++)
            {
                if (mesh.Intersections[i].FilledCount == 2)
                {
                    startIntersection = i;
                    break;
                }
            }
            if (startIntersection == -1)
                return false;

            int visitedCount = 0;
            int currentIntersection = startIntersection;
            int prevEdge = -1;

            while (visitedCount < filledEdgeCount)
            {
                int nextEdge = -1;
                foreach (int edgeIdx in mesh.Intersections[currentIntersection].Edges)
                {
                    if (edgeIdx != prevEdge && mesh.Edges[edgeIdx].State == EdgeState.Filled)
                    {
                        nextEdge = edgeIdx;
                        break;
                    }
                }
                if (nextEdge == -1)
                    break;

                visitedCount++;
                Edge edge = mesh.Edges[nextEdge];
                currentIntersection = edge.Intersections[0] == currentIntersection ? edge.Intersections[1] : edge.Intersections[0];
                prevEdge = nextEdge;
                if (currentIntersection == startIntersection)
                    break;
            }

            return visitedCount == filledEdgeCount && currentIntersection == startIntersection;
        }
    }
}
