#nullable enable
using System;
using System.Collections.Generic;

namespace TwistNTurn
{
    public static class PuzzleHelper
    {
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

        public static void GetDefaultSize(MeshType type, out int width, out int height)
        {
            switch (type)
            {
                case MeshType.Hexagonal3:
                    width = 3;
                    height = 3;
                    break;
                default:
                    width = 10;
                    height = 10;
                    break;
            }
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
                if (!int.TryParse(wSpan, out width) || !int.TryParse(hSpan, out height))
                    return false;
                return width > 0 && height > 0;
            }
            else
            {
                if (!int.TryParse(span, out width))
                    return false;
                height = width;
                return width > 0;
            }
        }

        public static MeshType MeshTypeFromString(string typeName)
        {
            if (string.IsNullOrEmpty(typeName))
                return MeshType.Square;
            return typeName.Replace(" ", "").Replace("-", "").Replace("_", "").ToLowerInvariant() switch
            {
                "square" or "sq" => MeshType.Square,
                "squaresymmetrical" or "sqs" => MeshType.SquareSymmetrical,
                "hexagon3" or "hexagonal3" or "hex3" => MeshType.Hexagonal3,
                _ => MeshType.Square
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

            var visitedEdges = new HashSet<int>();
            int currentIntersection = startIntersection;
            int prevEdge = -1;

            while (true)
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
                if (nextEdge == -1 || visitedEdges.Contains(nextEdge))
                {
                    if (nextEdge != -1) visitedEdges.Add(nextEdge);
                    break;
                }
                visitedEdges.Add(nextEdge);
                Edge edge = mesh.Edges[nextEdge];
                int nextInter = edge.Intersections[0] == currentIntersection ? edge.Intersections[1] : edge.Intersections[0];
                currentIntersection = nextInter;
                prevEdge = nextEdge;
            }

            return visitedEdges.Count == filledEdgeCount && currentIntersection == startIntersection;
        }
    }
}
