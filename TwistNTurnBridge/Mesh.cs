using System;
using System.Collections.Generic;
using System.Text;
#if !BRIDGE
using System.IO;
#endif
using System.Diagnostics;
using System.Linq;
using TwistNTurnBridge;

namespace TwistNTurn
{
    public enum EdgeState
    {
        Empty,
        Filled,
        Excluded,
    }

    public enum TriState
    {
        Unknown,
        Same,
        Opposite,
    }

    public enum IntersType
    {
        Unknown,
        Twist,
        Turn,
    }

    public enum EdgePairRestriction
    {
        None,
        NotBoth,
        NotNeither,
    }

    public enum SolverMethod
    {
        Iterative,
        Recursive,
    }

    public enum SolveState
    {
        NoSolutions,
        MultipleSolutions,
        Solved,
    }

    public class Edge
    {
        public Edge()
        {
            Cells = new List<int>();
            Intersections = new int[2];
        }
        public Edge(Edge other)
        {
            State = other.State;
            Intersections = (int[])other.Intersections.Clone();
            Cells = new List<int>(other.Cells);
            Color = other.Color;
            EdgeSet = other.EdgeSet;
        }

        public List<int> Cells;
        public int[] Intersections;
        public EdgeState State;
        public int Color;
        public int EdgeSet;

        public Edge Clone()
        {
            return new Edge(this);
        }
    }

    public class Intersection
    {
        public Intersection()
        {
            Edges = new List<int>();
            Opposite = new List<int>();
            Cells = new List<int>();
        }
        public Intersection(Intersection other)
        {
            FilledCount = other.FilledCount;
            ExcludedCount = other.ExcludedCount;
            Edges = new List<int>(other.Edges);
            Opposite = new List<int>(other.Opposite);
            Cells = new List<int>(other.Cells);
            Type = other.Type;
            X = other.X;
            Y = other.Y;
            EdgeSet = other.EdgeSet;
            // Can't clone edge set entries across, have to 'fixup' them later.
            //EdgeSetEntry = other.EdgeSetEntry;
        }
        public List<int> Edges;
        public List<int> Opposite;
        public List<int> Cells;
        public int FilledCount;
        public int ExcludedCount;
        public float X;
        public float Y;
        public int EdgeSet;
        public ChainNode EdgeSetEntry;
        public IntersType Type;

        public Intersection Clone()
        {
            return new Intersection(this);
        }
    }

    public class Cell
    {
        public Cell()
        {
            Edges = new List<int>();
            Intersections = new List<int>();
        }
        public Cell(Cell other)
        {
            Edges = new List<int>(other.Edges);
            Intersections = new List<int>(other.Intersections);
            Color = other.Color;
            FilledCount = other.FilledCount;
            ExcludedCount = other.ExcludedCount;
        }
        public List<int> Edges;
        public List<int> Intersections;
        public int FilledCount;
        public int ExcludedCount;
        public int Color;

        public Cell Clone()
        {
            return new Cell(this);
        }
    }

    public enum MeshType
    {
        Square,
        Hexagonal3,
        SquareSymmetrical,
    }

    #region ApproxPointStorage class to help with constructing grids.

    public class ApproxPointStorage
    {

        public ApproxPointStorage(float eps)
        {
            comparer = new EpsComparer(eps);
            lookup = new Dictionary<Point, int>(comparer);
        }
        EpsComparer comparer;

        Dictionary<Point, int> lookup;

        public int Add(float x, float y, int newIndex)
        {
            Point p = new Point(x, y);
            if (lookup.ContainsKey(p))
                return lookup[p];
            else
                lookup.Add(p, newIndex);
            return newIndex;
        }

        class Point
        {
            public Point(float x, float y)
            {
                X = x;
                Y = y;
            }
            public float X;
            public float Y;
        }

        class EpsComparer : IEqualityComparer<Point>
        {
            public EpsComparer(float eps)
            {
                this.eps = eps;
            }
            float eps;

            private float Round(float coord)
            {
                coord += (float) (Math.PI * Math.E);
                return coord - (coord % eps);
            }

            public bool Equals(Point x, Point y)
            {
                return Math.Abs(Round(x.X) - Round(y.X)) < eps/2 && Math.Abs(Round(x.Y) - Round(y.Y)) < eps/2;
            }

            public int GetHashCode(Point obj)
            {
                return (Round(obj.X).GetHashCode()*33) ^ Round(obj.Y).GetHashCode();
            }
        }
    }
    #endregion

    public delegate void MeshChangeUpdateEventHandler(object sender, MeshChangeUpdateEventArgs args);

    public class MeshChangeUpdateEventArgs
    {
        public MeshChangeUpdateEventArgs(Mesh current, List<IAction> justDone, bool success)
        {
            this.CurrentMesh = current;
            this.JustDone = justDone;
            this.SuccessfulAttempt = success;
        }
        public MeshChangeUpdateEventArgs(Mesh current, List<IAction> justDone, bool success, bool starting)
        {
            this.CurrentMesh = current;
            this.JustDone = justDone;
            this.SuccessfulAttempt = success;
            this.Starting = starting;
        }

        public bool SuccessfulAttempt;

        public List<IAction> JustDone;

        public Mesh CurrentMesh;

        public bool Starting;
    }

    public class Mesh
    {
        public string DebugView()
        {
            int bottomLeft = -1;
            for (int i = 1; i < intersections.Count; i++)
            {
                if (intersections[i].Edges.Count == 2)
                {
                    bottomLeft = i;
                    break;
                }
            }
            int height = bottomLeft + 1;
            int width = intersections.Count / height;
            char[] display = new char[4 * height * (width+1)];
            for (int i = 0; i < display.Length; i++)
                display[i] = ' ';
            for (int i = 0; i < intersections.Count; i++)
            {
                int x = i / height;
                int y = i % height;
                int index = 2 * x + (2 * (width * 2+1) * y);
                switch (intersections[i].Type)
                {
                    case IntersType.Unknown:
                        display[index] = '+';
                        break;
                    case IntersType.Turn:
                        display[index] = '*';
                        break;
                    case IntersType.Twist:
                        display[index] = '@';
                        break;
                }
            }
            for (int i = 0; i < height; i++)
            {
                display[2 * i * (2 * width + 1) + (2 * width + 1) - 1] = '\n';
                display[2 * i * (2 * width + 1) + (2 * width + 1) - 2] = '\r';
                display[(2 * i+1) * (2 * width + 1) + (2 * width + 1) - 1] = '\n';
                display[(2 * i+1) * (2 * width + 1) + (2 * width + 1) - 2] = '\r';
            }
            for (int i = 0; i < edges.Count; i++)
            {
                int span = 2 * height - 1;
                if (i >= span * (width - 1))
                {
                    int y = i - span * (width - 1);
                    display[2 * y * (2 * width + 1) + (2 * width + 1) + 2 * width - 2] = CharForState(edges[i].State, true);
                }
                else
                {
                    int x = i / span;
                    int y = i % span;
                    display[y * (2 * width + 1) + 2 * x + (y % 2 == 0 ? 1 : 0)] = CharForState(edges[i].State, y % 2 != 0);
                }
            }
            return new string(display);
        }

        private char CharForState(EdgeState edgeState, bool vert)
        {
            switch (edgeState)
            {
                case EdgeState.Empty:
                    return ' ';
                case EdgeState.Excluded:
                    return 'x';
                case EdgeState.Filled:
                    return vert ?'|' : '-';
            }
            return ' ';
        }

        public Mesh(Mesh other)
        {
            edges = new List<Edge>();
            foreach (Edge edge in other.edges)
                edges.Add(edge.Clone());
            intersections = new List<Intersection>();
            foreach (Intersection inters in other.intersections)
                intersections.Add(inters.Clone());
            cells = new List<Cell>();
            foreach (Cell cell in other.cells)
                cells.Add(cell.Clone());
            SolverMethod = other.SolverMethod;
            IterativeSolverDepth = other.iterativeSolverDepth;
            IterativeRecMaxDepth = other.iterativeRecMaxDepth;
            meshType = other.meshType;
            considerMultipleLoops = other.considerMultipleLoops;
            this.coloringCheats = other.coloringCheats;
            this.colorSets = new List<List<int>>();
            foreach (List<int> otherColorSet in other.colorSets)
                colorSets.Add(new List<int>(otherColorSet));
            this.cellColorSets = new List<List<int>>();
            foreach (List<int> otherCellColorSet in other.cellColorSets)
                cellColorSets.Add(new List<int>(otherCellColorSet));
            this.contaminateFullSolver = other.contaminateFullSolver;
            this.edgeSets = new List<List<int>>();
            foreach (List<int> otherEdgeSet in other.edgeSets)
                edgeSets.Add(new List<int>(otherEdgeSet));
            this.generateLengthFraction = other.generateLengthFraction;
            this.satisifiedIntersCount = other.satisifiedIntersCount;
            this.useColoring = other.useColoring;
            this.useCellPairs = other.useCellPairs;
            this.useCellPairsTopLevel = other.useCellPairsTopLevel;
            this.useEdgeRestricts = other.useEdgeRestricts;
        }

        #region Mesh Construction using known prototype.
        public Mesh(int width, int height, MeshType type)
        {
            meshType = type;
            edges = new List<Edge>();
            intersections = new List<Intersection>();
            cells = new List<Cell>();
            if (type == MeshType.Square || type == MeshType.SquareSymmetrical)
            {
                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                    {
                        Intersection newInters = new Intersection();
                        intersections.Add(newInters);
                        newInters.X = i;
                        newInters.Y = j;
                        if (i > 0 && j < height)
                            newInters.Cells.Add(j + (i - 1) * height);
                        if (i < width && j < height)
                            newInters.Cells.Add(j + i * height);
                        if (j > 0 && i < width)
                            newInters.Cells.Add((j - 1) + i * height);
                        if (i > 0 && j > 0)
                            newInters.Cells.Add((j - 1) + (i - 1) * height);
                        if (i < width && j < height)
                        {
                            Cell newCell = new Cell();
                            cells.Add(newCell);
                            newCell.Intersections.Add(j + i * (height + 1));
                            newCell.Intersections.Add(j + 1 + i * (height + 1));
                            newCell.Intersections.Add(j + 1 + (i + 1) * (height + 1));
                            newCell.Intersections.Add(j + (i + 1) * (height + 1));
                        }
                    }
                }
                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                    {
                        if (i < width)
                        {
                            int index = edges.Count;
                            Edge newEdge = new Edge();
                            edges.Add(newEdge);
                            newEdge.Intersections[0] = j + i * (height + 1);
                            newEdge.Intersections[1] = j + (i + 1) * (height + 1);
                            foreach (int interIndex in newEdge.Intersections)
                            {
                                intersections[interIndex].Edges.Add(index);
                            }
                            if (j > 0)
                                newEdge.Cells.Add(j - 1 + i * height);
                            if (j < height)
                                newEdge.Cells.Add(j + i * height);
                            foreach (int cellIndex in newEdge.Cells)
                            {
                                cells[cellIndex].Edges.Add(index);
                            }
                        }
                        if (j < height)
                        {
                            int index = edges.Count;
                            Edge newEdge = new Edge();
                            edges.Add(newEdge);
                            newEdge.Intersections[0] = j + i * (height + 1);
                            newEdge.Intersections[1] = j + 1 + i * (height + 1);
                            foreach (int interIndex in newEdge.Intersections)
                            {
                                intersections[interIndex].Edges.Add(index);
                            }
                            if (i > 0)
                                newEdge.Cells.Add(j + (i - 1) * height);
                            if (i < width)
                                newEdge.Cells.Add(j + i * height);
                            foreach (int cellIndex in newEdge.Cells)
                            {
                                cells[cellIndex].Edges.Add(index);
                            }
                        }
                    }
                }
                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                    {
                        if (i < width)
                        {
                            if (i == 0)
                            {
                                intersections[j + i * (height + 1)].Opposite.Add(-1);
                            }
                            else
                            {
                                intersections[j + i * (height + 1)].Opposite.Add(GetEdgeJoining(j + i * (height + 1), j + (i - 1) * (height + 1)));
                            }
                            if (i < width - 1)
                            {
                                intersections[j + (i + 1) * (height + 1)].Opposite.Add(GetEdgeJoining(j + (i + 1) * (height + 1), j + (i + 2) * (height + 1)));
                            }
                            else
                            {
                                intersections[j + (i + 1) * (height + 1)].Opposite.Add(-1);
                            }

                        }
                        if (j < height)
                        {
                            if (j == 0)
                            {
                                intersections[j + i * (height + 1)].Opposite.Add(-1);
                            }
                            else
                            {
                                intersections[j + i * (height + 1)].Opposite.Add(GetEdgeJoining(j + i * (height + 1), j - 1 + (i) * (height + 1)));
                            }
                            if (j < height - 1)
                            {
                                intersections[j + 1 + i * (height + 1)].Opposite.Add(GetEdgeJoining(j + 1 + (i) * (height + 1), j + 2 + (i) * (height + 1)));
                            }
                            else
                            {
                                intersections[j + 1 + i * (height + 1)].Opposite.Add(-1);
                            }
                        }
                    }
                }
            }
            else if (type == MeshType.Hexagonal3)
            {
                ApproxPointStorage storage = new ApproxPointStorage(0.001f);
                float hexagonWidth = 1 + 2 * (float)Math.Cos(Math.PI / 3);
                float hexagonHeight = 2 * (float)Math.Cos(Math.PI / 6);
                float scalingFactor = 1.5f;
                float hexagonWidthBit = (float)Math.Cos(Math.PI / 3) * scalingFactor;
                hexagonHeight *= scalingFactor;
                hexagonWidth *= scalingFactor;

                for (int i = 0; i < width; i++)
                {
                    for (int j = 0; j < height; j++)
                    {
                        List<int> intersects = new List<int>();
                        float topX = hexagonHeight * i + scalingFactor * i + (j % 2 == 1 ? hexagonHeight / 2 + scalingFactor / 2 : 0) + hexagonHeight + scalingFactor;
                        float topY = (float)Math.Sqrt(3) / 2 * (hexagonHeight + scalingFactor) * j + hexagonWidth;
                        AddIntersection(storage, intersects, topX - hexagonHeight / 2 - scalingFactor, topY + scalingFactor / 2);
                        AddIntersection(storage, intersects, topX - hexagonHeight / 2, topY + scalingFactor / 2);
                        AddIntersection(storage, intersects, topX - hexagonHeight / 2, topY - scalingFactor / 2);
                        AddIntersection(storage, intersects, topX - hexagonHeight / 2 - scalingFactor, topY - scalingFactor / 2);
                        AddPolyBoundry(intersects);
                        intersects.Clear();
                        AddIntersection(storage, intersects, topX + hexagonHeight / 2 + scalingFactor, topY + scalingFactor / 2);
                        AddIntersection(storage, intersects, topX + hexagonHeight / 2, topY + scalingFactor / 2);
                        AddIntersection(storage, intersects, topX + hexagonHeight / 2, topY - scalingFactor / 2);
                        AddIntersection(storage, intersects, topX + hexagonHeight / 2 + scalingFactor, topY - scalingFactor / 2);
                        AddPolyBoundry(intersects);
                        intersects.Clear();
                        AddIntersection(storage, intersects, topX - hexagonHeight / 2, topY - scalingFactor / 2);
                        AddIntersection(storage, intersects, topX - hexagonHeight / 2 - hexagonWidthBit, topY - scalingFactor / 2 - hexagonHeight / 2);
                        AddIntersection(storage, intersects, topX - hexagonWidthBit, topY - scalingFactor / 2 - hexagonHeight / 2 - hexagonWidthBit);
                        AddIntersection(storage, intersects, topX, topY - scalingFactor / 2 - hexagonWidthBit);
                        AddPolyBoundry(intersects);
                        intersects.Clear();
                        AddIntersection(storage, intersects, topX, topY - scalingFactor / 2 - hexagonWidthBit);
                        AddIntersection(storage, intersects, topX + hexagonWidthBit, topY - scalingFactor / 2 - hexagonHeight / 2 - hexagonWidthBit);
                        AddIntersection(storage, intersects, topX + hexagonHeight / 2 + hexagonWidthBit, topY - scalingFactor / 2 - hexagonHeight / 2);
                        AddIntersection(storage, intersects, topX + hexagonHeight / 2, topY - scalingFactor / 2);
                        AddPolyBoundry(intersects);
                        intersects.Clear();
                        AddIntersection(storage, intersects, topX - hexagonHeight / 2, topY + scalingFactor / 2);
                        AddIntersection(storage, intersects, topX - hexagonHeight / 2 - hexagonWidthBit, topY + scalingFactor / 2 + hexagonHeight / 2);
                        AddIntersection(storage, intersects, topX - hexagonWidthBit, topY + scalingFactor / 2 + hexagonHeight / 2 + hexagonWidthBit);
                        AddIntersection(storage, intersects, topX, topY + scalingFactor / 2 + hexagonWidthBit);
                        AddPolyBoundry(intersects);
                        intersects.Clear();
                        AddIntersection(storage, intersects, topX, topY + scalingFactor / 2 + hexagonWidthBit);
                        AddIntersection(storage, intersects, topX + hexagonWidthBit, topY + scalingFactor / 2 + hexagonHeight / 2 + hexagonWidthBit);
                        AddIntersection(storage, intersects, topX + hexagonHeight / 2 + hexagonWidthBit, topY + scalingFactor / 2 + hexagonHeight / 2);
                        AddIntersection(storage, intersects, topX + hexagonHeight / 2, topY + scalingFactor / 2);
                        AddPolyBoundry(intersects);
                        intersects.Clear();
                    }
                }
                CreateOpposites();
                CreateCells();
            }
            FullClear();
        }
        #endregion

        #region Mesh construction helpers.

        private void AddPolyBoundry(List<int> intersects)
        {
            for (int k = 0; k < intersects.Count; k++)
            {
                int start = intersects[k];
                int end = intersects[(k + 1) % intersects.Count];
                try
                {
                    GetEdgeJoining(start, end);
                }
                catch
                {
                    AddEdge(start, end);
                }
            }
        }

        private void AddIntersection(ApproxPointStorage storage, List<int> intersects, float topX, float topY)
        {
            int intersectNum = storage.Add(topX, topY, intersections.Count);
            intersects.Add(intersectNum);
            if (intersectNum == intersections.Count)
            {
                Intersection inters = new Intersection();
                intersections.Add(inters);
                inters.X = topX;
                inters.Y = topY;
            }
        }

        private void AddIntersectionOnGrid(int[,] grid, List<int> intersects, int topX, int topY)
        {
            if (grid[topX, topY] == -1)
            {
                grid[topX, topY] = intersections.Count;
                intersects.Add(intersections.Count);
                Intersection inters = new Intersection();
                intersections.Add(inters);
                inters.X = topX;
                inters.Y = topY;
            }
            else
                intersects.Add(grid[topX, topY]);
        }

        public void CreateCells()
        {
            for (int i = 0; i < intersections.Count; i++)
            {
                List<List<int>> cellCornersSets = CreateCells(i);
                foreach (List<int> cellCorners in cellCornersSets)
                {
                    if (cellCorners.Count > 0)
                    {
                        int cellIndex = cells.Count;
                        Cell cell = new Cell();
                        cells.Add(cell);
                        cell.Intersections.AddRange(cellCorners);
                        for (int j = 0; j < cellCorners.Count; j++)
                        {
                            int cellCorner = cellCorners[j];
                            Intersection inters = intersections[cellCorner];
                            inters.Cells.Add(cellIndex);
                            int eIndex = GetEdgeJoining(cellCorner, cellCorners[(j + 1) % cellCorners.Count]);
                            Edge e = edges[eIndex];
                            e.Cells.Add(cellIndex);
                            cell.Edges.Add(eIndex);
                        }
                    }
                }
            }
        }

        public void CreateOpposites()
        {
            for (int i = 0; i < intersections.Count; i++)
            {
                for (int j = 0; j < intersections[i].Edges.Count; j++)
                {
                    intersections[i].Opposite.Add(ReverseEdge(intersections[i].Edges[j], i));
                }
            }
        }

        private int ReverseEdge(int e, int i)
        {
            int otherI = GetOtherInters(e, i);
            float startX = intersections[i].X;
            float startY = intersections[i].Y;
            float endX = intersections[otherI].X;
            float endY = intersections[otherI].Y;
            double idealAngle = Math.Atan2(startY - endY, startX - endX);
            // In hexagon/square/triangle the 'ideal candidate' might be 30 degrees off.
            // So this gives a tolerance of 36 degrees to ensure it isn't missed, but reduce probability
            // of a scenario where we incorrectly choose an edge as opposite.
            // In practice both square and hexagon/square/triangle only have scenarios of 90 degrees being
            // the closest potential incorect candidate when there should be no candidate.
            double bestAngle = Math.PI / 5;
            int bestCandidate = -1;

            for (int j = 0; j < intersections[i].Edges.Count; j++)
            {
                if (intersections[i].Edges[j] == e) continue;
                int candidate = GetOtherInters(intersections[i].Edges[j], i);
                float candidateX = intersections[candidate].X;
                float candidateY = intersections[candidate].Y;
                double candidateAngle = Math.Atan2(candidateY - startY, candidateX - startX);
                double angleDiff = Math.Abs(idealAngle - candidateAngle);

                if (angleDiff >= Math.PI) angleDiff = 2*Math.PI - angleDiff;
                if (angleDiff < bestAngle)
                {
                    bestAngle = angleDiff;
                    bestCandidate = intersections[i].Edges[j];
                }
            }
            return bestCandidate;
            /*
            // TODO - support other than square.
            // TODO - optimize this to not result in quadratic performance in number of intersections.
            int startX = (int)intersections[i].X;
            int startY = (int)intersections[i].Y;
            int endX = (int)intersections[otherI].X;
            int endY = (int)intersections[otherI].Y;
            if (startY == endY)
            {
                int xdiff = endX - startX;
                int xtarget = startX - xdiff;
                int result = -1;
                for (int j = 0; j < intersections.Count; j++)
                {
                    if ((int) intersections[j].X == xtarget && (int) intersections[j].Y == startY)
                    {
                        result = j;
                        break;
                    }
                }
                return result == -1 ? -1 : GetEdgeJoining(i, result);
            }
            else
            {
                int ydiff = endY - startY;
                int ytarget = startY - ydiff;
                int result = -1;
                for (int j = 0; j < intersections.Count; j++)
                {
                    if ((int)intersections[j].X == startX && (int)intersections[j].Y == ytarget)
                    {
                        result = j;
                        break;
                    }
                }
                return result == -1 ? -1 : GetEdgeJoining(i, result);
            }
            */
        }

        private List<List<int>> CreateCells(int i)
        {
            List<List<int>> result = new List<List<int>>();
            Intersection inters = intersections[i];
            foreach (int eIndex in inters.Edges)
            {
                List<int> corners = new List<int>();
                int last = i;
                int next = GetOtherInters(eIndex, last);
                Intersection nextInter = intersections[next];
                Intersection lastInter = inters;
                double lastAngle = Math.Atan2(nextInter.Y - inters.Y, nextInter.X - inters.X);
                corners.Add(i);
                while (next != i)
                {
                    corners.Add(next);
                    int realLast = last;
                    last = next;
                    lastInter = nextInter;
                    next = -1;
                    // switch lastAngle into our space.
                    lastAngle = Math.PI + lastAngle;
                    if (lastAngle > Math.PI)
                        lastAngle -= Math.PI * 2;
                    double nextAngle = lastAngle + Math.PI * 2;
                    foreach (int edgeIndex in nextInter.Edges)
                    {
                        int possibleNext = GetOtherInters(edgeIndex, last);
                        if (possibleNext == realLast)
                            continue;
                        Intersection possibleNextInter = intersections[possibleNext];
                        double possibleNextAngle = Math.Atan2(possibleNextInter.Y - lastInter.Y, possibleNextInter.X - lastInter.X);
                        if (possibleNextAngle > lastAngle && possibleNextAngle < nextAngle)
                        {
                            nextAngle = possibleNextAngle;
                            next = possibleNext;
                        }
                        else
                        {
                            possibleNextAngle += Math.PI * 2;
                            if (possibleNextAngle > lastAngle && possibleNextAngle < nextAngle)
                            {
                                nextAngle = possibleNextAngle;
                                next = possibleNext;
                            }
                        }
                    }
                    if (nextAngle > Math.PI)
                        nextAngle -= Math.PI * 2;
                    lastAngle = nextAngle;
                    nextInter = intersections[next];
                }
                int min = int.MaxValue;
                foreach (int corner in corners)
                {
                    if (corner < min)
                        min = corner;
                }
                if (min == i)
                {
                    double total = 0.0;
                    // Verify anticlockwise.
                    for (int cIndex = 0; cIndex < corners.Count; cIndex++)
                    {
                        int start = corners[cIndex];
                        int mid = corners[(cIndex + 1) % corners.Count];
                        int end = corners[(cIndex + 2) % corners.Count];
                        Intersection startInter = intersections[start];
                        Intersection midInter = intersections[mid];
                        Intersection endInter = intersections[end];
                        float dy1 = midInter.Y - startInter.Y;
                        float dx1 = midInter.X - startInter.X;
                        float dy2 = endInter.Y - midInter.Y;
                        float dx2 = endInter.X - midInter.X;
                        double cross = dy1 * dx2 - dx1 * dy2;
                        double dot = dx1 * dx2 + dy1 * dy2;
                        double angle = Math.Acos(dot / Math.Sqrt(dx1 * dx1 + dy1 * dy1) / Math.Sqrt(dx2 * dx2 + dy2 * dy2));
                        if (cross < 0)
                            angle = -angle;

                        total += angle;
                    }
                    if (total > 0)
                        result.Add(corners);
                }
            }
            return result;
        }

        private int GetOtherInters(int eIndex, int last)
        {
            Edge e = edges[eIndex];
            if (e.Intersections[0] == last)
                return e.Intersections[1];
            else
                return e.Intersections[0];
        }

        public void AddEdge(int start, int end)
        {
            int index = edges.Count;
            Edge newEdge = new Edge();
            edges.Add(newEdge);
            newEdge.Intersections[0] = start;
            newEdge.Intersections[1] = end;
            foreach (int interIndex in newEdge.Intersections)
            {
                intersections[interIndex].Edges.Add(index);
            }

        }

        #endregion

        public event MeshChangeUpdateEventHandler MeshChangeUpdate;

        public MeshType MeshType
        {
            get
            {
                return meshType;
            }
        }
        private MeshType meshType;

        public List<Edge> Edges
        {
            get
            {
                return edges;
            }
        }
        List<Edge> edges;

        public List<Intersection> Intersections
        {
            get
            {
                return intersections;
            }
        }
        List<Intersection> intersections;

        public List<Cell> Cells
        {
            get
            {
                return cells;
            }
        }
        List<Cell> cells;

        private int[,] edgeDistances
        {
            get
            {
                if (edgeDistancesCache == null)
                {
                    edgeDistancesCache = new int[edges.Count, edges.Count];
                    for (int i = 0; i < edges.Count; i++)
                    {
                        for (int j = 0; j < edges.Count; j++)
                        {
                            if (i != j)
                                edgeDistancesCache[i, j] = int.MaxValue;
                        }
                    }
                    bool progress = true;
                    while (progress)
                    {
                        progress = false;
                        for (int next = 0; next < edges.Count; next++)
                        {
                            Edge e = edges[next];
                            List<int> oneAway = new List<int>();
                            foreach (int inters in e.Intersections)
                            {
                                Intersection i = intersections[inters];
                                foreach (int other in i.Edges)
                                {
                                    if (other != next)
                                        oneAway.Add(other);
                                }
                            }
                            for (int i = 0; i < edges.Count; i++)
                            {
                                int toMe = edgeDistancesCache[next, i];
                                if (toMe != int.MaxValue)
                                {
                                    foreach (int other in oneAway)
                                    {
                                        if (edgeDistancesCache[other, i] > toMe + 1)
                                        {
                                            edgeDistancesCache[other, i] = edgeDistancesCache[i, other] = toMe + 1;
                                            progress = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                return edgeDistancesCache;
            }
        }
        int[,] edgeDistancesCache;

        public void GetEdgeExtent(Edge e, out float startX, out float startY, out float endX, out float endY)
        {
            Intersection start = intersections[e.Intersections[0]];
            Intersection end = intersections[e.Intersections[1]];
            startX = start.X;
            startY = start.Y;
            endX = end.X;
            endY = end.Y;
        }

        public SolverMethod SolverMethod
        {
            get
            {
                return solverMethod;
            }
            set
            {
                solverMethod = value;
            }
        }
        private SolverMethod solverMethod;

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

        public Mesh FinalSolution
        {
            get
            {
                return finalSolution;
            }
        }
        Mesh finalSolution;

        public int[] FinalDepthPatern
        {
            get
            {
                return finalDepthPatern;
            }
        }
        int[] finalDepthPatern;

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

        public int GetEdgeJoining(int inters1, int inters2)
        {
            Intersection rInters1 = intersections[inters1];
            Intersection rInters2 = intersections[inters2];
            for (int i = 0; i < rInters1.Edges.Count; i++)
            {
                int edgeIndex = rInters1.Edges[i];
                for (int j = 0; j < rInters2.Edges.Count; j++)
                {
                    int otherEdgeIndex = rInters2.Edges[j];
                    if (edgeIndex == otherEdgeIndex)
                    {
                        return edgeIndex;
                    }
                }
            }
            throw new ArgumentException("Specified intersections aren't joined.");
        }

        private int GetOppositeEdge(int edge, int inters)
        {
            Intersection i = intersections[inters];
            for (int j = 0; j < i.Edges.Count; j++)
            {
                if (i.Edges[j] == edge)
                {
                    return i.Opposite[j];
                }
            }
            throw new InvalidOperationException();
        }


        public bool JoinColor(int a, int b, bool invert, List<int[]> changes, ref bool wasteOftime)
        {
            wasteOftime = false;
            Edge e = edges[a];
            Edge other = edges[b];
            if (e.Color == 0 && other.Color == 0)
            {
                colorSets.Add(new List<int>());
                changes.Add(new int[] { colorSets.Count });
                e.Color = colorSets.Count;
                colorSets[e.Color - 1].Add(a);
                changes.Add(new int[] { a, 0, e.Color });
                other.Color = invert ? -e.Color : e.Color;
                colorSets[e.Color - 1].Add(b);
                changes.Add(new int[] { b, 0, other.Color });
            }
            else if (e.Color == 0)
            {
                e.Color = invert ? -other.Color : other.Color;
                colorSets[Math.Abs(e.Color) - 1].Add(a);
                changes.Add(new int[] { a, 0, e.Color });
            }
            else if (other.Color == 0)
            {
                other.Color = invert ? -e.Color : e.Color;
                colorSets[Math.Abs(other.Color) - 1].Add(b);
                changes.Add(new int[] { b, 0, other.Color });
            }
            else
            {
                if (invert && e.Color == other.Color)
                    return false;
                if (!invert && e.Color == -other.Color)
                    return false;
                if (e.Color == other.Color)
                {
                    wasteOftime = true;
                    return true;
                }
                if (e.Color == -other.Color)
                {
                    wasteOftime = true;
                    return true;
                }
                int before = Math.Abs(e.Color) > Math.Abs(other.Color) ? e.Color : other.Color;
                int after = Math.Abs(e.Color) > Math.Abs(other.Color) ? other.Color : e.Color;
                if (invert)
                    after = -after;
                List<int> colorSetBefore = colorSets[Math.Abs(before) - 1];
                List<int> colorSetAfter = colorSets[Math.Abs(after) - 1];
                for (int i = colorSetBefore.Count - 1; i >= 0; i--)
                {
                    int edge = colorSetBefore[i];
                    Edge curEdge = edges[edge];
                    if (curEdge.Color == before)
                    {
                        curEdge.Color = after;
                        colorSetAfter.Add(edge);
                        changes.Add(new int[] { edge, before, after });
                    }
                    else
                    {
                        curEdge.Color = -after;
                        colorSetAfter.Add(edge);
                        changes.Add(new int[] { edge, -before, -after });
                    }
                }
                colorSets[Math.Abs(before) - 1].Clear();
            }
            return true;
        }


        public bool JoinCellColor(int a, int b, bool invert, List<int[]> changes, ref bool wasteOftime)
        {
            //Debug.WriteLine(string.Format("CellJoin {0} to {1} invert={2}", a, b, invert), "Actions");
            wasteOftime = false;
            Cell c = cells[a];
            if (b == -1)
            {
                if (cellColorSets.Count == 0)
                {
                    cellColorSets.Add(new List<int>());
                    changes.Add(new int[] { cellColorSets.Count });
                }
                if (!invert && c.Color == 1 || invert && c.Color == -1)
                {
                    wasteOftime = true;
                    return true;
                }
                if (Math.Abs(c.Color) == 1)
                {
                    UnjoinCellColor(changes);
                    return false;
                }
                if (c.Color == 0)
                {
                    c.Color = invert ? -1 : 1;
                    cellColorSets[0].Add(a);
                    changes.Add(new int[] { a, 0, c.Color });
                }
                else
                {
                    int before = c.Color;
                    int after = 1;
                    if (invert)
                        after = -after;
                    List<int> colorSetBefore = cellColorSets[Math.Abs(before) - 1];
                    List<int> colorSetAfter = cellColorSets[0];
                    for (int i = colorSetBefore.Count - 1; i >= 0; i--)
                    {
                        int cell = colorSetBefore[i];
                        Cell curCell = cells[cell];
                        if (curCell.Color == before)
                        {
                            curCell.Color = after;
                            colorSetAfter.Add(cell);
                            changes.Add(new int[] { cell, before, after });
                        }
                        else
                        {
                            curCell.Color = -after;
                            colorSetAfter.Add(cell);
                            changes.Add(new int[] { cell, -before, -after });
                        }
                    }
                    cellColorSets[Math.Abs(before) - 1].Clear();
                }
                return true;
            }
            Cell other = cells[b];
            if (c.Color == 0 && other.Color == 0)
            {
                cellColorSets.Add(new List<int>());
                // First color set is reserved for the non-existant cell.
                changes.Add(new int[] { cellColorSets.Count });
                if (cellColorSets.Count == 1)
                {
                    cellColorSets.Add(new List<int>());
                    changes.Add(new int[] { cellColorSets.Count });
                }
                c.Color = cellColorSets.Count;
                cellColorSets[c.Color - 1].Add(a);
                changes.Add(new int[] { a, 0, c.Color });
                other.Color = invert ? -c.Color : c.Color;
                cellColorSets[c.Color - 1].Add(b);
                changes.Add(new int[] { b, 0, other.Color });
            }
            else if (c.Color == 0)
            {
                c.Color = invert ? -other.Color : other.Color;
                cellColorSets[Math.Abs(c.Color) - 1].Add(a);
                changes.Add(new int[] { a, 0, c.Color });
            }
            else if (other.Color == 0)
            {
                other.Color = invert ? -c.Color : c.Color;
                cellColorSets[Math.Abs(other.Color) - 1].Add(b);
                changes.Add(new int[] { b, 0, other.Color });
            }
            else
            {
                if (invert && c.Color == other.Color)
                    return false;
                if (!invert && c.Color == -other.Color)
                    return false;
                if (c.Color == other.Color)
                {
                    wasteOftime = true;
                    return true;
                }
                if (c.Color == -other.Color)
                {
                    wasteOftime = true;
                    return true;
                }
                int before = Math.Abs(c.Color) > Math.Abs(other.Color) ? c.Color : other.Color;
                int after = Math.Abs(c.Color) > Math.Abs(other.Color) ? other.Color : c.Color;
                if (invert)
                    after = -after;
                List<int> colorSetBefore = cellColorSets[Math.Abs(before) - 1];
                List<int> colorSetAfter = cellColorSets[Math.Abs(after) - 1];
                for (int i = colorSetBefore.Count - 1; i >= 0; i--)
                {
                    int cell = colorSetBefore[i];
                    Cell curCell = cells[cell];
                    if (curCell.Color == before)
                    {
                        curCell.Color = after;
                        colorSetAfter.Add(cell);
                        changes.Add(new int[] { cell, before, after });
                    }
                    else
                    {
                        curCell.Color = -after;
                        colorSetAfter.Add(cell);
                        changes.Add(new int[] { cell, -before, -after });
                    }
                }
                cellColorSets[Math.Abs(before) - 1].Clear();
            }
            return true;
        }

        public double GenerateLengthFraction
        {
            get
            {
                return generateLengthFraction;
            }
            set
            {
                generateLengthFraction = value;
            }
        }
        private double generateLengthFraction = 0.5;

        public double GenerateBoringFraction
        {
            get
            {
                return generateBoringFraction;
            }
            set
            {
                generateBoringFraction = value;
            }
        }
        private double generateBoringFraction = 0.01;

        public bool CreatingLoop = false;

        public void Generate()
        {
            AbortPrune = false;
            bool done = false;
            List<IAction> backup = new List<IAction>();
            Random rnd = new Random();
            int targetCount = (int)Math.Floor(intersections.Count * generateLengthFraction);
            long countSum = 0;
            int tries = 0;
            while (!done)
            {
                try
                {
                    CreatingLoop = true;
                    done = true;
                    int loopTries = 0;
                    bool loopToSmall = true;
                    while (loopToSmall)
                    {
                        loopToSmall = false;
                        FullClear();
                        backup.Clear();
                        int start = rnd.Next(intersections.Count);
                        Intersection inters = intersections[start];
                        int edgeIntersIndex = rnd.Next(inters.Edges.Count);
                        int edgeIndex = inters.Edges[edgeIntersIndex];
                        Perform(edgeIndex, EdgeState.Filled, backup, 0);
                        bool success = CreateLoop(rnd, start, start, edgeIntersIndex);
                        int count = 0;
                        foreach (Edge edge in edges)
                            if (edge.State == EdgeState.Filled)
                                count++;
                        countSum += count;
                        tries++;
                        if (count < targetCount || RateBoringness() > GenerateBoringFraction)
                        {
                            if (loopTries > 100 && count >= countSum / tries)
                            {
                                if (TryExpandLoop(rnd, targetCount - count))
                                    break;
                            }
                            loopToSmall = true;
                            loopTries++;
                            if (loopTries > 1000)
                                break;
                        }
                    }
                    for (int i = 0; i < 8; i++)
                    {
                        UpdateTwists();
                        /*Mesh clone = new Mesh(this);
                        clone.Clear();
                        SolveState state = clone.TrySolve();
                        Mesh clone2 = new Mesh(this);
                        if (state == SolveState.Solved && MakeLoopBetter())
                            throw new Exception(clone2.DebugView());*/
                        
                        if (!MakeLoopBetter())
                            break;
                        if (i < 7)
                            ClearTwists();
                        else
                        {
                            ClearTwists();
                            UpdateTwists();
                        }
                    }
                }
                finally
                {
                    CreatingLoop = false;
                }
                Clear();
                try
                {
                    PruneCounts();
                }
                catch (Exception e)
                {
                    if (e.Message == "Can't solve it anyway")
                        done = false;
                    else
                        throw;
                }
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
            List<int> edgesSeen = new List<int>();
            List<int> intersSeen = new List<int>();
            List<IntersType> intersTypes = new List<IntersType>();
            List<bool> turned = new List<bool>();
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
            List<int> candidates = new List<int>();
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
            Dictionary<int, bool> candidateCheck = new Dictionary<int, bool>();
            foreach (int i in candidates)
                candidateCheck.Add(i, true);
            bool changed = false;
            foreach (int edge in candidates)
            {
                changed |= TryPop(edge,  candidateCheck);
            }
            return changed;
        }


        private bool TryPop(int edgeToPop, Dictionary<int, bool> candidateCheck)
        {
            Edge eToPop = edges[edgeToPop];
            foreach (int c in eToPop.Cells)
            {
                Cell cell = cells[c];
                if (!IsPopable(cell, candidateCheck))
                    continue;
                int oldCount = cell.FilledCount;
                foreach (int edge in cell.Edges)
                {
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
                                foreach (int edgeIndex in inter.Edges)
                                {
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

        private bool IsPopable(Cell cell, Dictionary<int, bool> candidateCheck)
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

        private bool AllCandidates(Cell cell, Dictionary<int, bool> candidateCheck)
        {
            foreach (int edge in cell.Edges)
            {
                if (edges[edge].State == EdgeState.Filled && !candidateCheck.ContainsKey(edge))
                    return false;
            }
            return true;
        }        

        bool[] boringEdges;

        private double RateBoringness()
        {
            if (boringEdges == null)
                boringEdges = new bool[edges.Count];
            else
                Array.Clear(boringEdges, 0, boringEdges.Length);
            int boringCount = 0;
            int total = edges.Count;
            for (int index = 0; index < edges.Count; index++)
            {
                Edge e = edges[index];
                if (e.State != EdgeState.Empty)
                    continue;
                bool found = false;
                foreach (int interIndex in e.Intersections)
                {
                    Intersection inter = intersections[interIndex];
                    foreach (int edgeIndex in inter.Edges)
                    {
                        if (edges[edgeIndex].State != EdgeState.Empty)
                        {
                            found = true;
                            break;
                        }
                    }
                    if (found)
                        break;
                }
                if (!found)
                {
                    boringCount++;
                    boringEdges[index] = true;
                }
            }
            return (double)boringCount / (double)total;
        }

        private bool TryExpandLoop(Random rnd, int targetAdditional)
        {
            int lastTarget = 0;
            while (lastTarget != targetAdditional)
            {
                lastTarget = targetAdditional;
                List<int> cellsToExpand = new List<int>();
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
                for (int i = 0; i < cellsToExpand.Count; i++)
                {
                    int first = rnd.Next(cellsToExpand.Count);
                    int second = rnd.Next(cellsToExpand.Count);
                    int tmp = cellsToExpand[first];
                    cellsToExpand[first] = cellsToExpand[second];
                    cellsToExpand[second] = tmp;
                }
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
            foreach (int interIndex in cell.Intersections)
            {
                Intersection inter = intersections[interIndex];
                foreach (int edge in inter.Edges)
                {
                    if (boringEdges[edge])
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
            foreach (int interIndex in cell.Intersections)
            {
                Intersection inter = intersections[interIndex];
                if (inter.FilledCount != 2)
                    continue;
                bool found = false;
                foreach (int edgeIndex in inter.Edges)
                {
                    Edge edge = edges[edgeIndex];
                    if (edge.State == EdgeState.Filled)
                    {
                        foreach (int cellIndex in edge.Cells)
                        {
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
            foreach (int edge in cell.Edges)
            {
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
                            foreach (int edgeIndex in inter.Edges)
                            {
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
            Intersection curInters = null;
            int curIntersIndex = -1;
            foreach (int intersPos in prevEdge.Intersections)
            {
                Intersection other = intersections[intersPos];
                if (other != prevInters)
                {
                    curIntersIndex = intersPos;
                    curInters = other;
                    break;
                }
            }
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
                for (int i = 0; i < 2 * availEdges.Count; i++)
                {
                    int a = rnd.Next(availEdges.Count);
                    int b = rnd.Next(availEdges.Count);
                    if (a == b)
                        continue;
                    int tmp = availEdges[a];
                    availEdges[a] = availEdges[b];
                    availEdges[b] = tmp;
                }
            }
            List<IAction> backup = new List<IAction>();
            foreach (int availEdge in availEdges)
            {
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
                foreach (int edgeIndex in inters.Edges)
                {
                    Edge edge = edges[edgeIndex];
                    if (edge.State == EdgeState.Empty)
                    {
                        foreach (int otherInters in edge.Intersections)
                        {
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

        private IntersType DetermineType(int inter)
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
        private EdgePairRestriction[,] edgePairRestrictions;


        public event EventHandler PrunedCountProgress;

        public bool AbortPrune = false;

        private bool pruning = false;

        private void PruneCounts()
        {
            SolveState state = TrySolve();
            if (state != SolveState.Solved)
                throw new Exception("Can't solve it anyway");
            try
            {
                pruning = true;
                finalSolution = solutionsFound[0];
                finalDepthPatern = solutionDepthPatern[0];
                bool[] tried = new bool[intersections.Count];
                int[] trials = new int[intersections.Count];
                for (int i = 0; i < trials.Length; i++)
                    trials[i] = i;
                Random rnd = new Random();
                for (int i = 0; i < trials.Length * 2; i++)
                {
                    int a = rnd.Next(trials.Length);
                    int b = rnd.Next(trials.Length);
                    int tmp = trials[a];
                    trials[a] = trials[b];
                    trials[b] = tmp;
                }
                int width = 0;
                int height = 0;
                if (meshType == MeshType.SquareSymmetrical)
                {
                    bool found = true;
                    int iSearch = 0;
                    while (found)
                    {
                        found = true;
                        try
                        {
                            GetEdgeJoining(iSearch, iSearch + 1);
                        }
                        catch
                        {
                            found = false;
                        }
                        iSearch++;
                    }
                    height = iSearch - 1;
                    width = Cells.Count / height;

                }
                foreach (int trial in trials)
                {
                    if (AbortPrune)
                        break;
                   // if (meshType != MeshType.SquareSymmetrical)
                   // {
                        Intersection inter = intersections[trial];
                        IntersType oldVal = inter.Type;
                        inter.Type = IntersType.Unknown;
                        if (TrySolve() != SolveState.Solved)
                        {
                            inter.Type = oldVal;
                        }
                        else
                        {
                            finalSolution = solutionsFound[0];
                            finalDepthPatern = solutionDepthPatern[0];
                        }
                        if (PrunedCountProgress != null)
                            PrunedCountProgress(this, EventArgs.Empty);
                    //}
                    //else
                    //{
                    //    Cell cell = cells[trial];
                    //    if (tried[trial])
                    //    {
                    //        continue;
                    //    }
                    //    int oldVal = cell.TargetCount;
                    //    RemoveTarget(cell);
                    //    int y = trial / width;
                    //    int x = trial % width;
                    //    int otherY = height - y - 1;
                    //    int otherX = width - x - 1;
                    //    int otherTrial = otherY * width + otherX;
                    //    if (otherTrial == trial)
                    //    {
                    //        if (TrySolve() != SolveState.Solved)
                    //        {
                    //            AddTarget(cell, oldVal);
                    //        }
                    //        else
                    //        {
                    //            finalSolution = solutionsFound[0];
                    //            finalDepthPatern = solutionDepthPatern[0];
                    //        }
                    //        tried[trial] = true;
                    //        if (PrunedCountProgress != null)
                    //            PrunedCountProgress(this, EventArgs.Empty);
                    //    }
                    //    else
                    //    {
                    //        // While we have just removed a cell, there is no way trial and other trial can be adjacent
                    //        // Not in square symmetrical at least - so we don't need to recalculate the protected list.
                    //        if (prot.Contains(otherTrial))
                    //        {
                    //            AddTarget(cell, oldVal);
                    //            if (PrunedCountProgress != null)
                    //                PrunedCountProgress(this, EventArgs.Empty);
                    //            continue;
                    //        }
                    //        Cell otherCell = cells[otherTrial];
                    //        int otherOldVal = otherCell.TargetCount;
                    //        RemoveTarget(otherCell);
                    //        if (TrySolve() != SolveState.Solved)
                    //        {
                    //            AddTarget(cell, oldVal);
                    //            AddTarget(otherCell, otherOldVal);
                    //        }
                    //        else
                    //        {
                    //            finalSolution = solutionsFound[0];
                    //            finalDepthPatern = solutionDepthPatern[0];
                    //        }
                    //        tried[trial] = true;
                    //        tried[otherTrial] = true;
                    //        if (PrunedCountProgress != null)
                    //            PrunedCountProgress(this, EventArgs.Empty);
                    //        if (PrunedCountProgress != null)
                    //            PrunedCountProgress(this, EventArgs.Empty);
                    //    }
                    //}
                }
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
#if !BRIDGE
            catch (System.Threading.ThreadAbortException ex)
            {
                // We're screwed, no point trying to rollback changes.
                edgeChanges.Clear();
                return SolveState.NoSolutions;
            }
#endif
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
#if !BRIDGE
            catch (System.Threading.ThreadAbortException ex)
            {
                // We're screwed, no point trying to rollback changes.
                realChanges.Clear();
                return SolveState.NoSolutions;
            }
#endif
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
#if !BRIDGE
            catch (System.Threading.ThreadAbortException ex)
            {
                // We're screwed anyway, no point trying to undo the changes.
                realChanges.Clear();
                return SolveState.NoSolutions;
            }
#endif
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
        private Dictionary<int, int> oldFirstOrder;
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

        private bool IterativeTrySolveInternal(List<IAction> realChanges, ref bool changed, int iterativeRecDepth, IAction locusAction)
        {
            if (iterativeSolverDepth < 0)
                return true;
            bool iterTopLevel = true;
            List<IAction> trials = new List<IAction>();
            if (iterativeRecDepth != this.iterativeRecDepth)
            {
                iterTopLevel = false;
                GatherLocals(realChanges, trials, locusAction);
            }
            else
            {
                GatherAll(trials);
            }
            List<IAction> edgeChanges1 = new List<IAction>();
            List<IAction> edgeChanges2 = new List<IAction>();
            IAction lastTrial = null;
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
            foreach (IAction action in realChanges)
            {
                if (action is SetAction)
                {
                    SetAction setAction = (SetAction)action;
                    foreach (int edge in setAction.GetAffectedEdges())
                    {
                        GatherNearEdge(trials, edge, locusAction, true);
                    }
                }
                else if (action is ColorJoinAction)
                {
                    ColorJoinAction cjAction = (ColorJoinAction)action;
                    foreach (int edge in cjAction.GetAffectedEdges())
                    {
                        GatherNearEdge(trials, edge, locusAction, true);
                    }
                }
                else if (action is CellColorJoinAction)
                {
                    CellColorJoinAction ccjAction = (CellColorJoinAction)action;
                    foreach (int cell in ccjAction.GetAffectedCells())
                    {
                        foreach (int edge in cells[cell].Edges)
                        {
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
            foreach (int i in intersections[e.Intersections[0]].Edges)
            {
                if (edges[i].State == EdgeState.Empty && GetEdgeDistance(i, locusAction) <= maxDist)
                    smarts.Add(new SetAction(this, i, EdgeState.Filled));
            }
            foreach (int i in intersections[e.Intersections[1]].Edges)
            {
                if (edges[i].State == EdgeState.Empty && GetEdgeDistance(i, locusAction) <= maxDist)
                    smarts.Add(new SetAction(this, i, EdgeState.Filled));
            }
            if (useCellColoring && useCellColoringTrials)
            {
                foreach (int i in e.Cells)
                {
                    if (Math.Abs(cells[i].Color) != 1 && GetCellDistance(i, locusAction) <= maxDist)
                        smarts.Add(new CellColorJoinAction(this, i, -1, true));
                }
            }
            if (followCells)
            {
                foreach (int i in cells[e.Cells[0]].Edges)
                {
                    if (edges[i].State == EdgeState.Empty && GetEdgeDistance(i, locusAction) <= maxDist)
                        smarts.Add(new SetAction(this, i, EdgeState.Filled));
                }
                if (e.Cells.Count > 1)
                {
                    foreach (int i in cells[e.Cells[1]].Edges)
                    {
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

        EdgeState[] edgesSeen;
        TriState[,] edgePairsSeen;
        List<KeyValuePair<int, int>> edgePairsToClean = new List<KeyValuePair<int, int>>();
        TriState[,] cellPairsSeen;
        List<KeyValuePair<int, int>> cellPairsToClean = new List<KeyValuePair<int, int>>();
        EdgePairRestriction[,] edgeRestrictsSeen;
        List<KeyValuePair<int, int>> edgeRestrictsToClean = new List<KeyValuePair<int, int>>();
        bool[] cellsSeen;
        bool[] cellColorEdgeColorsSeen;
        bool[] intersectsSeen;
        bool[] interactsSeen;
        bool[] interactsSeen2;
        bool[] interactsSeen3;

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
        private bool Perform(IAction action, List<IAction> backup)
        {
            return Perform(action, backup, int.MaxValue);
        }
        private bool Perform(IAction action, List<IAction> backup, int maxDepth)
        {
            return Perform(action, backup, maxDepth, null);
        }

        List<IAction>[] moves;
        List<int>[] toConsiderEdges;
        List<int>[] toConsiderEdgeSets;
        List<int>[] toConsiderEdgeColors;
        List<int>[] toConsiderIntersects;
        List<int>[] toConsiderCellColors;

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

        private bool Perform(IAction action, List<IAction> backup, int maxDepth, List<int> initialIntersectsAffected)
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

            List<int> colorSetsSeen = new List<int>();
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
                foreach (IAction move in moves[curDepth])
                {
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
                foreach (IAction move in moves[curDepth])
                {
                    if (move is SetAction)
                    {
                        SetAction setAction = (SetAction)move;
                        AddEdgetToConsider(action, toConsiderEdges, curDepth, setAction.EdgeIndex);
                        foreach (int edge in setAction.GetAffectedEdges())
                        {
                            AddEdgetToConsider(action, toConsiderEdgeSets, curDepth, edge);
                        }
                    }
                    else if (move is ColorJoinAction)
                    {
                        ColorJoinAction cjAction = (ColorJoinAction)move;
                        foreach (int edge in cjAction.GetAffectedEdges())
                        {
                            AddEdgetToConsider(action, toConsiderEdgeColors, curDepth, edge);
                        }
                    }
                    else if (move is CellColorJoinAction)
                    {
                        CellColorJoinAction ccjAction = (CellColorJoinAction)move;
                        foreach (int cell in ccjAction.GetAffectedCells())
                        {
                            AddCellToConsider(action, toConsiderCellColors, curDepth, cell);
                        }

                    }
                    else if (move is EdgeRestrictionAction)
                    {
                        EdgeRestrictionAction erAction = (EdgeRestrictionAction)move;
                        foreach (int edge in erAction.GetAffectedEdges())
                        {
                            AddEdgetToConsider(action, toConsiderEdgeColors, curDepth, edge);
                        }
                    }
                }
                toConsiderEdges[curDepth].Sort();
                toConsiderEdgeSets[curDepth].Sort();
                toConsiderEdgeColors[curDepth].Sort();
                toConsiderIntersects[curDepth].Sort();
                toConsiderCellColors[curDepth].Sort();

                if (!ConsiderEdges(moves, toConsiderEdges, curDepth, colorSetsSeen))
                    return false;
                if (!ConsiderIntersects(moves, toConsiderIntersects, curDepth))
                    return false;
                if (!ConsiderEdgeSets(moves, toConsiderEdgeSets, curDepth))
                    return false;
                if (!ConsiderEdgeColors(moves, toConsiderEdgeColors, colorSetsSeen, curDepth))
                    return false;
                if (!ConsiderCellColors(moves, toConsiderCellColors, curDepth))
                    return false;

            }
            return true;
        }

        private void ClearCellPairs()
        {
            foreach (KeyValuePair<int, int> kvp in cellPairsToClean)
            {
                cellPairsSeen[kvp.Key, kvp.Value] = TriState.Unknown;
                cellPairsSeen[kvp.Value, kvp.Key] = TriState.Unknown;
            }
            cellPairsToClean.Clear();
        }

        private void ClearEdgeRestricts()
        {
            foreach (KeyValuePair<int, int> kvp in edgeRestrictsToClean)
            {
                edgeRestrictsSeen[kvp.Key, kvp.Value] = EdgePairRestriction.None;
                edgeRestrictsSeen[kvp.Value, kvp.Key] = EdgePairRestriction.None;
            }
            edgeRestrictsToClean.Clear();
        }

        private void ClearEdgePairs()
        {
            foreach (KeyValuePair<int, int> kvp in edgePairsToClean)
            {
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
                foreach (int cellColorIndex in toConsiderCellColors[curDepth])
                {
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
                        if (!GatherCellColoringEdgeColoringMovesForCellColorChange(cell, moves, curDepth, cellColorIndex))
                            return false;
                    }
                }
            }
            return true;
        }

        private bool ConsiderEdgeColors(List<IAction>[] moves, List<int>[] toConsiderEdgeColors, List<int> colorSetsSeen, int curDepth)
        {
            if (UseColoring || UseEdgeRestricts)
            {
                int lastEdge = -1;
                foreach (int edgeAffectedIndex in toConsiderEdgeColors[curDepth])
                {
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
                        if (!GatherInteractForcedMoves(inters, moves, curDepth, edgesSeen, edgePairsSeen, edgeRestrictsSeen))
                            return false;
                    }
                    if (edge.Color != 0)
                    {
                        if (!GatherFollowColoringColorSetChanged(moves, curDepth, edge, edgesSeen, edgeAffectedIndex, colorSetsSeen))
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
                foreach (int edgeAffectedIndex in toConsiderEdgeSets[curDepth])
                {
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
            foreach (int intersectIndex in toConsiderIntersects[curDepth])
            {
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
                
                if ( UseColoring || UseEdgeRestricts)
                {
                    if (!interactsSeen2[intersectIndex])
                    {
                        interactsSeen2[intersectIndex] = true;
                        if (!GatherInteractForcedMoves(inters, moves, curDepth, edgesSeen, edgePairsSeen, edgeRestrictsSeen))
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

        private bool ConsiderEdges(List<IAction>[] moves, List<int>[] toConsiderEdges, int curDepth, List<int> colorSetsSeen)
        {
            int lastEdge = -1;
            foreach (int edgeAffectedIndex in toConsiderEdges[curDepth])
            {
                if (edgeAffectedIndex == lastEdge)
                    continue;
                lastEdge = edgeAffectedIndex;

                Edge edge = edges[edgeAffectedIndex];

                foreach (int intersIndex in edge.Intersections)
                {
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
                    foreach (int intersIndex in edge.Intersections)
                    {
                        if (interactsSeen2[intersIndex])
                            continue;
                        else
                            interactsSeen2[intersIndex] = true;
                        Intersection inters = intersections[intersIndex];
                        if (!GatherInteractForcedMoves(inters, moves, curDepth, edgesSeen, edgePairsSeen, edgeRestrictsSeen))
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
                        if (!GatherFollowColoring(moves, curDepth, edge, edgesSeen, edgeAffectedIndex, colorSetsSeen))
                            return false;
                    }
                }
                if (!GatherFollowEdgeRestrictions(moves, curDepth, edge, edgesSeen, edgeAffectedIndex))
                    return false;
            }
            return true;
        }

        private void AddEdgetToConsider(IAction sourceAction, List<int>[] toConsider, int curDepth, int other)
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

        private void AddCellToConsider(IAction sourceAction, List<int>[] toConsider, int curDepth, int cell)
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
            foreach (int intersIndex in edge.Intersections)
            {
                Intersection inters = intersections[intersIndex];
                if (inters.FilledCount != 1)
                    continue;
                foreach (int otherEdgeIndex in inters.Edges)
                {
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

        private bool GatherFollowColoringColorSetChanged(List<IAction>[] moves, int curDepth, Edge edge, EdgeState[] edgesSeen, int edgeIndex, List<int> colorSetsSeen)
        {
            List<int> colorSet;
            colorSet = colorSets[Math.Abs(edge.Color) - 1];
            // Stop caching that we have seen this if we were.
            colorSetsSeen.Remove(Math.Abs(edge.Color) - 1);
            EdgeState posState = EdgeState.Empty;
            foreach (int i in colorSet)
            {
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
        private bool GatherFollowColoring(List<IAction>[] moves, int curDepth, Edge edge, EdgeState[] edgesSeen, int edgeIndex, List<int> colorSetsSeen)
        {
            List<int> colorSet;
            colorSet = colorSets[Math.Abs(edge.Color) - 1];
            if (colorSetsSeen.Contains(Math.Abs(edge.Color) - 1))
                return true;
            colorSetsSeen.Add(Math.Abs(edge.Color) - 1);
            foreach (int i in colorSet)
            {
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
            if (!UseEdgeRestricts)
                return true;
            if (edge.State == EdgeState.Empty)
                return true;
            foreach (int inters in edge.Intersections)
            {
                Intersection inter = intersections[inters];
                foreach (int otherEdgeIndex in inter.Edges)
                {
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
            foreach (int cellIndex in edge.Cells)
            {
                Cell cell = cells[cellIndex];
                foreach (int otherEdgeIndex in cell.Edges)
                {
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
            if (edgePairRestrictions[edge1, edge2] != EdgePairRestriction.None)
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
                foreach (int otherEdgeIndex in inters.Edges)
                {
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
                foreach (int edge in cell.Edges)
                {
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
                Cell cell = null;
                int cellIndex = -2;
                Cell otherCell = null;
                int otherC = -2;
                foreach (int c in e.Cells)
                {
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
                if (otherC != -2)
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
                    List<int> cellColorSet1 = otherCellColor != 0 ? cellColorSets[Math.Abs(otherCellColor) - 1] : null;
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
                int hashcode = curNumber.GetHashCode();
                for (int i = 0; i < success.Length; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= success[i].GetHashCode();
                }
                return hashcode;
            }

            public override bool Equals(object obj)
            {
                SuccessLookup other = (SuccessLookup)obj;
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
            int[,] map;
            SuccessLookup key = new SuccessLookup{success=success, curNumber=curNumber};
            if (!successLookup.TryGetValue(key, out map))
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
            List<int[]> result = RetrieveActions(targets, baseLine, curNumber, numbering);
            if (result == null)
                return false;
            return ProcessRetrievedActions(moves, curDepth, edgesSeen, edgeNumber, edgePairsSeen, edgeRestrictsSeen, result);
      }

        private bool ProcessRetrievedActions(List<IAction>[] moves, int curDepth, EdgeState[] edgesSeen, int[] edgeNumber, TriState[,] edgePairsSeen, EdgePairRestriction[,] edgeRestrictsSeen, List<int[]> result)
        {
            foreach (int[] action in result)
            {
                if (action[0] == 0)
                {
                    int j = action[1];
                    bool reallyFilled = action[2] == 1;
                    Edge e = edges[edgeNumber[j]];
                    if (e.State == EdgeState.Empty)
                    {
                        if (!AddSetAction(edgeNumber[j], (reallyFilled ? EdgeState.Filled : EdgeState.Excluded), moves, curDepth + 1, edgesSeen))
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
                    EdgePairRestriction restr = action[3] == 1 ? EdgePairRestriction.NotBoth : EdgePairRestriction.NotNeither;
                    if (!AddEdgeRestrictAction(edgeNumber[n], edgeNumber[m], restr, moves, curDepth + 1, edgeRestrictsSeen, edgePairsSeen))
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
                    hashcode ^=  targets[i].Key.GetHashCode();
                    for (int j = 0; j < targets[i].Value.Count; j++)
                    {
                        hashcode += hashcode << 5;
                        hashcode ^= targets[i].Value[j].GetHashCode();
                    }
                }
                for (int i = 0; i < baseLine.Length; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= baseLine[i].GetHashCode();
                }
                for (int i = 0; i < numbering.Length; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= numbering[i].GetHashCode();
                }
                for (int i = 0; i < edgePatterns.Count; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= edgePatterns[i].GetHashCode();
                }
                for (int i = 0; i < edgeMasks.Count; i++)
                {
                    hashcode += hashcode << 5;
                    hashcode ^= edgeMasks[i].GetHashCode();
                }
                hashcode += hashcode << 5;
                hashcode ^= curNumber.GetHashCode();
                return hashcode;
            }
            public override bool Equals(object obj)
            {
                PatternLookup other = (PatternLookup)obj;
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

        Dictionary<PatternLookup, List<int[]>> patternLookup = new Dictionary<PatternLookup, List<int[]>>();
        List<uint> emptyList = new List<uint>();

        private List<int[]> RetrieveActions(List<KeyValuePair<uint, List<int>>> targets, int[] baseLine, int curNumber, int[] numbering)
        {
            List<int[]> result;
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
                    int[,] maps = null;
                    if (UseColoring || UseEdgeRestricts)
                        maps = GetMaps(success, curNumber);
                    if (UseColoring)
                    {
                        for (int i = 0; i < curNumber - 1; i++)
                        {
                            for (int j = i + 1; j < curNumber; j++)
                            {
                                int value = maps[i, j];
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
                                int value = maps[i, j];
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



        public bool LoadFromText(string[] lines)
        {
            // might have a single intersection, get rid of it...
            Intersections.Clear();
            int curLine = 1;
            if (lines[curLine] != "Intersections")
                return false;
            curLine++;
            int intersectionCount = int.Parse(lines[curLine]);
            curLine++;
            for (int i = 0; i < intersectionCount; i++)
            {
                string[] splits = lines[curLine].Split(' ');
                Intersection inters = new Intersection();
                inters.X = float.Parse(splits[0]);
                inters.Y = float.Parse(splits[1]);
                inters.Type = (IntersType)Enum.Parse(typeof(IntersType), splits[2], true);
                Intersections.Add(inters);
                curLine++;
            }
            if (lines[curLine] != "Edges")
                return false;
            curLine++;
            int edgeCount = int.Parse(lines[curLine]);
            curLine++;
            List<IAction> settings = new List<IAction>();
            for (int i = 0; i < edgeCount; i++)
            {
                string[] splits = lines[curLine].Split(' ');
                AddEdge(int.Parse(splits[0]), int.Parse(splits[1]));
                EdgeState state = (EdgeState)Enum.Parse(typeof(EdgeState), splits[2], true);
                if (state != EdgeState.Empty)
                {
                    settings.Add(new SetAction(this, GetEdgeJoining(int.Parse(splits[0]), int.Parse(splits[1])), state));
                }
                curLine++;
            }
            CreateOpposites();
            CreateCells();
            Clear();
            if (curLine < lines.Length && lines[curLine] == "EdgeColorSets")
            {
                curLine++;
                int edgeColorSetsCount = int.Parse(lines[curLine]);
                curLine++;
                for (int i = 0; i < edgeColorSetsCount; i++)
                {
                    curLine++;
                    int edgeColorSetCount = int.Parse(lines[curLine]);
                    curLine++;
                    int first = -1;
                    int firstColor = 0;
                    for (int j = 0; j < edgeColorSetCount; j++)
                    {
                        string[] bits = lines[curLine].Split(' ');
                        int edge = int.Parse(bits[0]);
                        int color = int.Parse(bits[1]);
                        if (first == -1)
                        {
                            first = edge;
                            firstColor = color;
                        }
                        else
                        {
                            settings.Add(new ColorJoinAction(this, first, edge, color == firstColor));
                        }
                        curLine++;
                    }
                }
            }
            if (curLine < lines.Length && lines[curLine] == "CellColorSets")
            {
                curLine++;
                int cellColorSetsCount = int.Parse(lines[curLine]);
                curLine++;
                for (int i = 0; i < cellColorSetsCount; i++)
                {
                    curLine++;
                    int cellColorSetCount = int.Parse(lines[curLine]);
                    curLine++;
                    int first = -1;
                    int firstColor = 0;
                    for (int j = 0; j < cellColorSetCount; j++)
                    {
                        string[] bits = lines[curLine].Split(' ');
                        int cell = int.Parse(bits[0]);
                        int color = int.Parse(bits[1]);
                        if (Math.Abs(color) != 1)
                        {
                            if (first == -1)
                            {
                                first = cell;
                                firstColor = color;
                            }
                            else
                            {
                                settings.Add(new CellColorJoinAction(this, first, cell, color == firstColor));
                            }
                        }
                        else
                        {
                            settings.Add(new CellColorJoinAction(this, cell, -1, color == 1));
                        }
                        curLine++;
                    }
                }
            }
            if (curLine < lines.Length && lines[curLine] == "EdgePairRestrictions")
            {
                curLine++;
                int edgePairRestrictionCount = int.Parse(lines[curLine]);
                curLine++;
                for (int j = 0; j < edgePairRestrictionCount; j++)
                {
                    string[] bits = lines[curLine].Split(' ');
                    int edge1 = int.Parse(bits[0]);
                    int edge2 = int.Parse(bits[1]);
                    EdgePairRestriction restriction = (EdgePairRestriction)Enum.Parse(typeof(EdgePairRestriction), bits[2], true);
                    settings.Add(new EdgeRestrictionAction(this, edge1, edge2, restriction));
                    curLine++;
                }
            }
            PerformListRegardless(settings);
            for (int i = 0; i < intersections.Count; i++)
            {
                if (intersections[i].Edges.Count != intersections[i].Opposite.Count)
                    throw new InvalidOperationException("WTH?");
            }
            return true;
        }

#if !BRIDGE
        public void Save(TextWriter writer)
        {
            writer.WriteLine(MeshType.ToString());
            writer.WriteLine("Intersections");
            writer.WriteLine(Intersections.Count);
            foreach (Intersection inters in Intersections)
            {
                writer.Write(inters.X);
                writer.Write(" ");
                writer.Write(inters.Y);
                writer.Write(" ");
                writer.Write(inters.Type);
                writer.WriteLine();
            }
            writer.WriteLine("Edges");
            writer.WriteLine(Edges.Count);
            foreach (Edge edge in Edges)
            {
                writer.Write(edge.Intersections[0]);
                writer.Write(" ");
                writer.Write(edge.Intersections[1]);
                writer.Write(" ");
                writer.WriteLine(edge.State.ToString());
            }
            writer.WriteLine("EdgeColorSets");
            writer.WriteLine(colorSets.Count);
            foreach (List<int> colorSet in colorSets)
            {
                writer.WriteLine("EdgeColorSet");
                writer.WriteLine(colorSet.Count);
                foreach (int edge in colorSet)
                {
                    writer.Write(edge);
                    writer.Write(" ");
                    writer.WriteLine(edges[edge].Color);
                }
            }
            writer.WriteLine("CellColorSets");
            writer.WriteLine(cellColorSets.Count);
            foreach (List<int> colorSet in cellColorSets)
            {
                writer.WriteLine("CellColorSet");
                writer.WriteLine(colorSet.Count);
                foreach (int cell in colorSet)
                {
                    writer.Write(cell);
                    writer.Write(" ");
                    writer.WriteLine(cells[cell].Color);
                }
            }
            writer.WriteLine("EdgePairRestrictions");
            int counter = 0;
            for (int i = 0; i < edges.Count; i++)
                for (int j = i + 1; j < edges.Count; j++)
                    if (edgePairRestrictions[i, j] != EdgePairRestriction.None)
                        counter++;
            writer.WriteLine(counter);
            for (int i = 0; i < edges.Count; i++)
                for (int j = i + 1; j < edges.Count; j++)
                    if (edgePairRestrictions[i, j] != EdgePairRestriction.None)
                    {
                        writer.Write(i);
                        writer.Write(" ");
                        writer.Write(j);
                        writer.Write(" ");
                        writer.WriteLine(edgePairRestrictions[i, j].ToString());
                    }
        }
#else
        public void Save(StringBuilder writer)
        {
            writer.AppendLine(MeshType.ToString());
            writer.AppendLine("Intersections");
            writer.AppendLine(Intersections.Count.ToString());
            foreach (Intersection inters in Intersections)
            {
                writer.Append(inters.X);
                writer.Append(" ");
                writer.Append(inters.Y);
                writer.Append(" ");
                writer.Append(inters.Type);
                writer.AppendLine();
            }
            writer.AppendLine("Edges");
            writer.AppendLine(Edges.Count.ToString());
            foreach (Edge edge in Edges)
            {
                writer.Append(edge.Intersections[0]);
                writer.Append(" ");
                writer.Append(edge.Intersections[1]);
                writer.Append(" ");
                writer.AppendLine(edge.State.ToString());
            }
            writer.AppendLine("EdgeColorSets");
            writer.AppendLine(colorSets.Count.ToString());
            foreach (List<int> colorSet in colorSets)
            {
                writer.AppendLine("EdgeColorSet");
                writer.AppendLine(colorSet.Count.ToString());
                foreach (int edge in colorSet)
                {
                    writer.Append(edge);
                    writer.Append(" ");
                    writer.AppendLine(edges[edge].Color.ToString());
                }
            }
            writer.AppendLine("CellColorSets");
            writer.AppendLine(cellColorSets.Count.ToString());
            foreach (List<int> colorSet in cellColorSets)
            {
                writer.AppendLine("CellColorSet");
                writer.AppendLine(colorSet.Count.ToString());
                foreach (int cell in colorSet)
                {
                    writer.Append(cell);
                    writer.Append(" ");
                    writer.AppendLine(cells[cell].Color.ToString());
                }
            }
            writer.AppendLine("EdgePairRestrictions");
            int counter = 0;
            for (int i = 0; i < edges.Count; i++)
            for (int j = i + 1; j < edges.Count; j++)
                if (edgePairRestrictions[i, j] != EdgePairRestriction.None)
                    counter++;
            writer.AppendLine(counter.ToString());
            for (int i = 0; i < edges.Count; i++)
            for (int j = i + 1; j < edges.Count; j++)
                if (edgePairRestrictions[i, j] != EdgePairRestriction.None)
                {
                    writer.Append(i);
                    writer.Append(" ");
                    writer.Append(j);
                    writer.Append(" ");
                    writer.AppendLine(edgePairRestrictions[i, j].ToString());
                }
        }
#endif


        internal bool PerformSetZero(int edgeIndex, EdgeState state, List<int[]> edgeSetChanges)
        {
            Edge edge = edges[edgeIndex];
            bool failed = RawEdgeSet(state, edge);
            if (state == EdgeState.Filled)
            {
                int edgeSet1 = GetEdgeSet(edge.Intersections[0], edgeIndex);
                int edgeSet2 = GetEdgeSet(edge.Intersections[1], edgeIndex);
                if (edgeSet1 == 0 && edgeSet2 == 0)
                {
                    edgeSets.Add(new List<int>());
                    edge.EdgeSet = edgeSets.Count;
                    edgeSets[edge.EdgeSet - 1].Add(edgeIndex);
                    edgeSetChanges.Add(new int[] { edgeIndex, edge.EdgeSet, 0 });
                }
                else if (edgeSet1 == 0 || edgeSet2 == 0)
                {
                    // one of them is positive
                    int edgeSet = Math.Max(edgeSet1, edgeSet2);
                    edgeSets[edgeSet - 1].Add(edgeIndex);
                    edge.EdgeSet = edgeSet;
                    edgeSetChanges.Add(new int[] { edgeIndex, edge.EdgeSet, 0 });
                }
                else if (edgeSet1 != edgeSet2)
                {
                    // Merge
                    int toKeep = Math.Min(edgeSet1, edgeSet2);
                    int toGo = Math.Max(edgeSet1, edgeSet2);
                    edgeSets[toKeep - 1].Add(edgeIndex);
                    edge.EdgeSet = toKeep;
                    edgeSetChanges.Add(new int[] { edgeIndex, edge.EdgeSet, 0 });
                    for (int i = edgeSets[toGo - 1].Count - 1; i >= 0; i--)
                    {
                        int otherEdge = edgeSets[toGo - 1][i];
                        edgeSets[toKeep - 1].Add(otherEdge);
                        Edge otherE = edges[otherEdge];
                        otherE.EdgeSet = toKeep;
                        edgeSetChanges.Add(new int[] { otherEdge, otherE.EdgeSet, toGo });
                    }
                    edgeSets[toGo - 1].Clear();
                }
                else
                {
                    int edgeSet = edgeSet1;
                    edgeSets[edgeSet - 1].Add(edgeIndex);
                    edge.EdgeSet = edgeSet;
                    edgeSetChanges.Add(new int[] { edgeIndex, edge.EdgeSet, 0 });
                    if (!failed && considerMultipleLoops)
                    {
                        // We're joining ... this might be bad.

                        // rule 1, can't close the loop if there are any intersections with odd number of filleds touching.
                        if (satisifiedIntersCount != intersections.Count)
                            return false;
                        // Rule 2, can't close the loop if there is more then one (non-empty) edge set.
                        int nonEmpty = 0;
                        foreach (List<int> otherEdgeSet in edgeSets)
                            if (otherEdgeSet.Count > 0)
                                nonEmpty++;
                        if (nonEmpty != 1)
                            return false;
                        // rule 3, can't close the loop if intersections aren't satifised.
                        // TODO: implemeht this more efficiently.
                        for (int i = 0; i < intersections.Count; i++)
                        {
                            if (intersections[i].Type != IntersType.Unknown)
                            {
                                if (intersections[i].Type != DetermineType(i))
                                    return false;
                            }
                        }
                        GenerateCheck();
                    }

                }
            }
            // up until here has to be atomic, all of it runs, or none of it does.  This is to ensure consistent state for rollback, later.
            if (failed)
                return false;
            return true;
        }

        private void GenerateCheck()
        {
            if (pruning && !earlyFail)
            {
                // We've found a solution, we can check that against the known solution to see if we've removed enough cells to force multiple solutions to exist.
                for (int i = 0; i < edges.Count; i++)
                {
                    if (edges[i].State == EdgeState.Filled)
                    {
                        if (finalSolution.edges[i].State != EdgeState.Filled)
                        {
                            // missmatch - early out somehow...
                            earlyFail = true;
                            break;
                        }
                    }
                }
            }
        }
        private bool earlyFail = false;

        private bool RawEdgeSet(EdgeState state, Edge edge)
        {
            if (edge.State != EdgeState.Empty)
                throw new Exception("Edge already set.");
            edge.State = state;
            bool failed = false;
            foreach (int cellIndex in edge.Cells)
            {
                Cell cell = cells[cellIndex];
                if (state == EdgeState.Filled)
                {
                    cell.FilledCount++;
                }
                else
                {
                    cell.ExcludedCount++;
                }
            }
            // TODO: consider failing immediately for straights on turns and turns on twists.
            foreach (int intersIndex in edge.Intersections)
            {
                Intersection inters = intersections[intersIndex];
                if (state == EdgeState.Filled)
                {
                    inters.FilledCount++;
                    if (inters.FilledCount > 2)
                        failed = true;
                    if (inters.FilledCount > 0 && inters.Edges.Count - inters.ExcludedCount < 2)
                        failed = true;
                    if (inters.FilledCount == 2)
                        satisifiedIntersCount++;
                    else if (inters.FilledCount < 4)
                        satisifiedIntersCount--;
                }
                else
                {
                    inters.ExcludedCount++;
                    if ((inters.FilledCount > 0 || inters.Type != IntersType.Unknown) && inters.Edges.Count - inters.ExcludedCount < 2)
                        failed = true;
                }
            }
            return failed;
        }

        private int GetEdgeSet(int intersection, int edgeToIgnore)
        {
            Intersection inter = intersections[intersection];
            foreach (int edgeIndex in inter.Edges)
            {
                if (edgeIndex == edgeToIgnore)
                    continue;
                Edge e = edges[edgeIndex];
                if (e.State == EdgeState.Filled)
                {
                    return e.EdgeSet;
                }
            }
            return 0;
        }

        internal void UnperformSetZero(int edgeIndex, EdgeState state, List<int[]> edgeSetChanges)
        {
            Edge edge = edges[edgeIndex];
            RawEdgeUnset(state, edge);
            for (int i = edgeSetChanges.Count - 1; i >= 0; i--)
            {
                int[] edgeSetChange = edgeSetChanges[i];
                int otherEdgeIndex = edgeSetChange[0];
                int newEdgeSet = edgeSetChange[1];
                int oldEdgeSet = edgeSetChange[2];

                List<int> edgeSet = edgeSets[newEdgeSet - 1];
                int edgeCheck = edgeSet[edgeSet.Count - 1];
                if (edgeCheck != otherEdgeIndex)
                    throw new InvalidOperationException("Attempting to undo out of order.");
                edgeSet.RemoveAt(edgeSet.Count - 1);
                if (oldEdgeSet != 0)
                {
                    edgeSet = edgeSets[oldEdgeSet - 1];
                    edgeSet.Add(otherEdgeIndex);
                }
                else
                {
                    // If it was a new edge, it may have been a new edge set, which needs to be removed.
                    if (newEdgeSet == edgeSets.Count && edgeSets[newEdgeSet - 1].Count == 0)
                        edgeSets.RemoveAt(newEdgeSet - 1);
                }
                Edge otherE = edges[otherEdgeIndex];
                otherE.EdgeSet = oldEdgeSet;
            }
        }

        private void RawEdgeUnset(EdgeState state, Edge edge)
        {
            edge.State = EdgeState.Empty;
            foreach (int cellIndex in edge.Cells)
            {
                Cell cell = cells[cellIndex];
                if (state == EdgeState.Filled)
                {
                    cell.FilledCount--;
                }
                else
                    cell.ExcludedCount--;
            }
            foreach (int intersIndex in edge.Intersections)
            {
                Intersection inters = intersections[intersIndex];
                if (state == EdgeState.Filled)
                {
                    inters.FilledCount--;
                    if (inters.FilledCount == 2 || inters.FilledCount == 0)
                        satisifiedIntersCount++;
                    else if (inters.FilledCount == 1)
                        satisifiedIntersCount--;
                }
                else
                    inters.ExcludedCount--;
            }
        }

        internal bool PerformUnsetZero(int edgeIndex, out EdgeState oldState, List<int[]> edgeSetChanges)
        {
            Edge edge = edges[edgeIndex];
            oldState = edge.State;
            RawEdgeUnset(oldState, edge);
            if (edge.EdgeSet != 0)
            {
                int oldEdgeSet = edge.EdgeSet;
                int index = edgeSets[edge.EdgeSet - 1].IndexOf(edgeIndex);
                edgeSetChanges.Add(new int[] { edgeIndex, 0, edge.EdgeSet, index });
                edgeSets[edge.EdgeSet - 1].RemoveAt(index);
                edge.EdgeSet = 0;
                int curInter = edge.Intersections[0];
                int target = edge.Intersections[1];
                List<int> foundEdges = new List<int>();
                Stack<int> todoEdges = new Stack<int>();
                Stack<int> todoPrevInters = new Stack<int>();
                todoEdges.Push(edgeIndex);
                todoPrevInters.Push(target);
                bool foundTarget = false;
                Dictionary<int, bool> foundInters = new Dictionary<int, bool>();
                while (todoEdges.Count > 0)
                {
                    int curEdge = todoEdges.Pop();
                    int lastInters = todoPrevInters.Pop();
                    int curInters = GetOtherInters(curEdge, lastInters);
                    if (curInters == target)
                    {
                        foundTarget = true;
                        break;
                    }
                    if (!foundInters.ContainsKey(curInters))
                    {
                        Intersection inters = intersections[curInters];
                        foreach (int nextEdge in inters.Edges)
                        {
                            if (nextEdge != curEdge)
                            {
                                Edge e = edges[nextEdge];
                                if (e.State == EdgeState.Filled)
                                {
                                    todoEdges.Push(nextEdge);
                                    todoPrevInters.Push(curInters);
                                }
                            }
                        }
                        foundInters[curInters] = true;
                    }
                }
                if (!foundTarget)
                {
                    edgeSets.Add(new List<int>());
                    int newEdgeSet = edgeSets.Count;
                    foreach (int otherEdgeIndex in foundEdges)
                    {
                        int otherIndex = edgeSets[oldEdgeSet - 1].IndexOf(otherEdgeIndex);
                        edgeSetChanges.Add(new int[] { otherEdgeIndex, newEdgeSet, oldEdgeSet, otherIndex });
                        edgeSets[oldEdgeSet - 1].RemoveAt(otherIndex);
                        Edge e = edges[otherEdgeIndex];
                        e.EdgeSet = newEdgeSet;
                        edgeSets[newEdgeSet - 1].Add(otherEdgeIndex);
                    }
                }
            }
            return true;
        }

        private int GetNextEdge(int curInter, int lastEdge)
        {
            Intersection inter = intersections[curInter];
            foreach (int edgeIndex in inter.Edges)
            {
                if (edgeIndex == lastEdge)
                    continue;
                Edge e = edges[edgeIndex];
                if (e.State == EdgeState.Filled)
                    return edgeIndex;
            }
            return -1;
        }

        internal void UnperformUnsetZero(int edgeIndex, EdgeState state, List<int[]> edgeSetChanges)
        {
            Edge edge = edges[edgeIndex];
            RawEdgeSet(state, edge);
            for (int i = edgeSetChanges.Count - 1; i >= 0; i--)
            {
                int[] edgeSetChange = edgeSetChanges[i];
                int otherEdgeIndex = edgeSetChange[0];
                int newEdgeSet = edgeSetChange[1];
                int oldEdgeSet = edgeSetChange[2];
                if (newEdgeSet != 0)
                {

                    List<int> edgeSet = edgeSets[newEdgeSet - 1];
                    int edgeCheck = edgeSet[edgeSet.Count - 1];
                    if (edgeCheck != otherEdgeIndex)
                        throw new InvalidOperationException("Attempting to undo out of order.");
                    edgeSet.RemoveAt(edgeSet.Count - 1);
                    if (oldEdgeSet != 0)
                    {
                        edgeSet = edgeSets[oldEdgeSet - 1];
                        edgeSet.Insert(edgeSetChange[3], otherEdgeIndex);
                    }
                    // If it was a new edge, it may have been a new edge set, which needs to be removed.
                    if (newEdgeSet == edgeSets.Count && edgeSets[newEdgeSet - 1].Count == 0)
                        edgeSets.RemoveAt(newEdgeSet - 1);
                }
                else
                {
                    List<int> edgeSet = edgeSets[oldEdgeSet - 1];
                    edgeSet.Insert(edgeSetChange[3], edgeIndex);
                }
                Edge otherE = edges[otherEdgeIndex];
                otherE.EdgeSet = oldEdgeSet;
            }

        }

        internal void UnjoinColor(List<int[]> colorSetChanges)
        {
            for (int i = colorSetChanges.Count - 1; i >= 0; i--)
            {
                if (colorSetChanges[i].Length == 1)
                {
                    if (colorSetChanges[i][0] == colorSets.Count && colorSets[colorSets.Count - 1].Count == 0)
                        colorSets.RemoveAt(colorSets.Count - 1);
                    else
                        throw new Exception("Color changes being performed out of order causing color leakage.");
                }
                else
                {
                    int edge = colorSetChanges[i][0];
                    int before = colorSetChanges[i][1];
                    int after = colorSetChanges[i][2];
                    Edge e = edges[edge];
                    if (e.Color != after)
                        throw new Exception("Colors being undone out of order.");
                    e.Color = before;
                    List<int> currentSet = colorSets[Math.Abs(after) - 1];
                    if (currentSet[currentSet.Count - 1] != edge)
                        throw new Exception("Colors being undone out of order.");
                    currentSet.RemoveAt(currentSet.Count - 1);
                    if (before != 0)
                    {
                        List<int> revertToSet = colorSets[Math.Abs(before) - 1];
                        revertToSet.Add(edge);
                    }
                }
            }
        }

        internal void UnjoinCellColor(List<int[]> colorSetChanges)
        {
            for (int i = colorSetChanges.Count - 1; i >= 0; i--)
            {
                if (colorSetChanges[i].Length == 1)
                {
                    //Debug.WriteLine(string.Format("Unjoin Cell dropping {0}", colorSetChanges[i][0]), "Actions");
                    if (colorSetChanges[i][0] == cellColorSets.Count && cellColorSets[cellColorSets.Count - 1].Count == 0)
                        cellColorSets.RemoveAt(cellColorSets.Count - 1);
                    else
                        throw new Exception("Color changes being performed out of order causing color leakage.");
                }
                else
                {
                    //Debug.WriteLine(string.Format("Unjoin Cell reverting {0} from {1} to {2}", colorSetChanges[i][0], colorSetChanges[i][2], colorSetChanges[i][1]), "Actions");
                    int cell = colorSetChanges[i][0];
                    int before = colorSetChanges[i][1];
                    int after = colorSetChanges[i][2];
                    Cell c = cells[cell];
                    if (c.Color != after)
                        throw new Exception("Colors being undone out of order.");
                    c.Color = before;
                    List<int> currentSet = cellColorSets[Math.Abs(after) - 1];
                    if (currentSet[currentSet.Count - 1] != cell)
                        throw new Exception("Colors being undone out of order.");
                    currentSet.RemoveAt(currentSet.Count - 1);
                    if (before != 0)
                    {
                        List<int> revertToSet = cellColorSets[Math.Abs(before) - 1];
                        revertToSet.Add(cell);
                    }
                }
            }
        }

        internal void GetSomething(List<IAction> changes)
        {
            IterativeTrySolveWithoutRollback(changes);
        }

        internal bool ClearCellColor(int cell1, out int index)
        {
            index = -1;
            Cell c = cells[cell1];
            if (c.Color == 0)
                return false;
            int colorSet = Math.Abs(c.Color) - 1;
            List<int> set = cellColorSets[colorSet];
            for (int i = 0; i < set.Count; i++)
            {
                if (set[i] == cell1)
                {
                    c.Color = 0;
                    index = i;
                    set.RemoveAt(i);
                    return true;
                }
            }
            return false;
        }

        internal void ResetCellColor(int cell1, int oldColor, int oldSetPos)
        {
            Cell c = cells[cell1];
            if (c.Color != 0)
                throw new InvalidOperationException("Unperforming actions out of order.");
            c.Color = oldColor;
            int colorSet = Math.Abs(c.Color) - 1;
            List<int> set = cellColorSets[colorSet];
            set.Insert(oldSetPos, cell1);
        }

        internal bool RestrictEdges(int edge1, int edge2, EdgePairRestriction state, ref bool wasteOfTime)
        {
            wasteOfTime = false;
            if (edgePairRestrictions[edge1, edge2] == state)
            {
                wasteOfTime = true;
                return true;
            }
            if (edgePairRestrictions[edge1, edge2] != EdgePairRestriction.None)
            {
                // Not matching is not a contradiction, technically we could 'upconvert' to a color join, but that is way too complicated for me right now.
                // TODO: upconvert to color join. WasteOfTime would then only possibly be false.
                wasteOfTime = true;
                return true;
            }
            edgePairRestrictions[edge1, edge2] = state;
            edgePairRestrictions[edge2, edge1] = state;
            return true;
        }

        internal void UnRestrictEdges(int edge1, int edge2)
        {
            Debug.Assert(this.edgePairRestrictions[edge1, edge2] != EdgePairRestriction.None);
            edgePairRestrictions[edge1, edge2] = EdgePairRestriction.None;
            edgePairRestrictions[edge2, edge1] = EdgePairRestriction.None;
        }

        Dictionary<int, List<KeyValuePair<int, int>>> edgeToPathSegmentMap = new Dictionary<int, List<KeyValuePair<int, int>>>();
        List<LoopPath> paths = new List<LoopPath>();
        int[,] pathStarts;

        internal bool TypeSatisfied(int i)
        {
            return intersections[i].Type == IntersType.Unknown || intersections[i].Type == DetermineType(i);
        }

        internal static Mesh FromPuzzleString(string puzzle)
        {
            string[] bits = puzzle.Split(':');
            if (bits.Length != 2)
                throw new ArgumentException("puzzle string invalid.");
            string[] dimensions = bits[0].Split('x');
            if (dimensions.Length != 2)
                throw new ArgumentException("puzzle string invalid.");
            int width;
            int height;
            if (!int.TryParse(dimensions[0], out width) || width < 2)
                throw new ArgumentException("puzzle string invalid.");
            if (!int.TryParse(dimensions[1], out height) || height < 2)
                throw new ArgumentException("puzzle string invalid.");
            Mesh result = new Mesh(width-1, height-1, MeshType.Square);
            int counter=0;
            foreach (char c in bits[1]) {
                int x = counter % width;
                int y = counter / width;
                int pos = y + x * height;
                if (c== '0') {
                    result.Intersections[pos].Type = IntersType.Twist;
                    counter++;
                }
                else if (c == '1')
                {
                    result.Intersections[pos].Type = IntersType.Turn;
                    counter++;
                }
                else
                {
                    counter += (int)(c - 'a' + 1);
                }
            }
            return result;
        }

        internal string ToPuzzleString()
        {
            StringBuilder result = new StringBuilder();

            int bottomLeft = -1;
            for (int i = 1; i < intersections.Count; i++)
            {
                if (intersections[i].Edges.Count == 2)
                {
                    bottomLeft = i;
                    break;
                }
            }
            int height = bottomLeft + 1;
            int width = intersections.Count / height;
            result.AppendFormat("{0}x{1}:", width, height);
            int counter = 0;
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int pos = y + x * height;
                    switch (intersections[pos].Type)
                    {
                        case IntersType.Twist:
                            {
                                if (counter > 0)
                                    result.Append((char)('a' + counter-1));
                                result.Append('0');
                                counter = 0;
                            }
                            break;
                        case IntersType.Turn:
                            {
                                if (counter > 0)
                                    result.Append((char)('a' + counter-1));
                                result.Append('1');
                                counter = 0;
                            }
                            break;
                        case IntersType.Unknown:
                            {
                                counter++;
                                if (counter > 23)
                                {
                                    result.Append((char)('a' + counter - 1));
                                    counter = 0;
                                }
                            }
                            break;
                    }
                }
            }
            return result.ToString();
        }

        internal string ToEdgeStateString()
        {
            return ToEdgeStateString(true);
        }

        internal string ToEdgeStateString(bool forceCross)
        {
            StringBuilder result = new StringBuilder();

            int bottomLeft = -1;
            for (int i = 1; i < intersections.Count; i++)
            {
                if (intersections[i].Edges.Count == 2)
                {
                    bottomLeft = i;
                    break;
                }
            }
            int height = bottomLeft + 1;
            int width = intersections.Count / height;
            int counter = 0;
            int counter2 = 0;
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int pos = y + x * height;
                    if (x + 1 < width)
                    {
                        int edge = GetEdgeJoining(pos, pos + height);
                        if (edges[edge].State == EdgeState.Filled)
                        {
                            AppendCounter(ref counter, result);
                            AppendCounter2(ref counter2, result);
                            result.Append('1');
                        }
                        else if (forceCross || edges[edge].State == EdgeState.Excluded)
                        {
                            AppendCounter2(ref counter2, result);
                            counter++;
                        }
                        else
                        {
                            AppendCounter(ref counter, result);
                            counter2++;
                        }
                    }
                    if (y + 1 < height)
                    {
                        int edge = GetEdgeJoining(pos, pos + 1);
                        if (edges[edge].State == EdgeState.Filled)
                        {
                            AppendCounter(ref counter, result);
                            AppendCounter2(ref counter2, result);
                            result.Append('1');
                        }
                        else if (forceCross || edges[edge].State == EdgeState.Excluded)
                        {
                            AppendCounter2(ref counter2, result);
                            counter++;
                        }
                        else
                        {
                            AppendCounter(ref counter, result);
                            counter2++;
                        }
                    }
                }
            }
            return result.ToString();
        }

        private void AppendCounter(ref int counter, StringBuilder result)
        {
            while (counter > 24)
            {
                result.Append((char)('a' + 24 - 1));
                counter -= 24;
            }
            if (counter > 0)
            {
                result.Append((char)('a' + counter - 1));
                counter = 0;
            }
        }
        private void AppendCounter2(ref int counter, StringBuilder result)
        {
            while (counter > 24)
            {
                result.Append((char)('A' + 24 - 1));
                counter -= 24;
            }
            if (counter > 0)
            {
                result.Append((char)('A' + counter - 1));
                counter = 0;
            }
        }

        internal void FromEdgeStateString(string input)
        {
            int bottomLeft = -1;
            for (int i = 1; i < intersections.Count; i++)
            {
                if (intersections[i].Edges.Count == 2)
                {
                    bottomLeft = i;
                    break;
                }
            }
            int height = bottomLeft + 1;
            int width = intersections.Count / height;
            int y = 0;
            int x = 0;
            bool horiz = true;
            Clear();
            foreach (char c in input)
            {
                int todo = 1;
                EdgeState toSet = EdgeState.Filled;
                if (char.IsLetter(c))
                {
                    if (char.IsLower(c))
                    {
                        toSet = EdgeState.Excluded;
                        todo = c - 'a';
                    }
                    else
                    {
                        toSet = EdgeState.Empty;
                        todo = c - 'A';
                    }
                }
                for (int i = 0; i < todo; i++)
                {
                    if (horiz && x + 1 >= width)
                    {
                        horiz = false;
                    }
                    int pos = y + x * height;
                    int edge;
                    if (horiz)
                    {
                        edge = GetEdgeJoining(pos, pos + height);
                    }
                    else
                    {
                        edge = GetEdgeJoining(pos, pos + 1);
                    }
                    if (toSet != EdgeState.Empty)
                        (new SetAction(this, edge, toSet)).Perform();
                    if (horiz && y + 1 < height)
                        horiz = false;
                    else if (x < width)
                    {
                        x++;
                        horiz = true;
                    }
                    else
                    {
                        y++;
                        x = 0;
                        horiz = true;
                    }
                }
            }
        }
    }

    // Path discovery - pattern check, adds a path.
    // Adding an edge which doesn't touch any other edges, adds a path.
    // Trial options, advanced merging logic can determine a path.  These paths can contain nested paths. (But if nested path is only proven on one trial we haven't proven a path...)
    // Adding an edge which touches an existing edge, extends a path. If the path is nested, the parent path entry needs updating.
    // Adding an exclusion which is in a path, excludes that path option.  
    // All edges of a path option are locked, exception nested paths.
    // If the last path option is excluded, all edges of the path are proven.
    // if two paths touch, all edges at that intersection which are not in either path are disproven.
    // if a path is added which touches two paths, perform validation.
    // Validation follows paths to determine the cycle.  
    // If a cycle exists then verifies that every path is in at least one path option, and that the longest path option is greater than the 'minimal path length'.
    // Additional validation, Verifies that no definitive paths occur as children of path options in such a way which proves a contradiction.
    // Specifically, any path which only occur as a child of one path, must not be in a different possibility branch to any other such path.
    // Alternatively, any path which occurs as a child of one path must be assumed to be the right path option, 
    // if this leads to the exclusion of a path from the overal path set, a contradiction has been found and validation has failed.

    public class LoopPath
    {
        public List<PathSegment> PathOptions;
        public int PathStart;
        public int PathEnd;
    }

    public class PathSegment
    {
        public LoopPath Parent;
        // Positive for an edge number, negative for a path index.
        public List<int> PathBits;
    }

    public class SetAction : IAction
    {

        public SetAction(Mesh mesh, int edge, EdgeState newState)
        {
            this.mesh = mesh;
            this.edge = edge;
            this.newState = newState;
        }

        private Mesh mesh;

        public int EdgeIndex
        {
            get
            {
                return edge;
            }
        }
        private int edge;

        public EdgeState EdgeState
        {
            get
            {
                return newState;
            }
        }
        private EdgeState newState;

        public List<int> GetAffectedEdges()
        {
            List<int> res = new List<int>();
            foreach (int[] change in edgeSetChanges)
            {
                // Ignore changes to the number of sets.
                if (change.Length > 1 && change[0] != edge)
                    res.Add(change[0]);
            }
            res.Add(edge);
            return res;
        }
        private List<int[]> edgeSetChanges;

        public bool Successful
        {
            get
            {
                return successful;
            }
        }
        private bool successful;

#region IAction Members

        public string Name
        {
            get { return "Make Edge " + edge.ToString() + " " + newState.ToString(); }
        }

        public bool Perform()
        {
            successful = true;
            edgeSetChanges = new List<int[]>();
            if (!mesh.PerformSetZero(edge, newState, edgeSetChanges))
            {
                successful = false;
            }
            return true;
        }

        public void Unperform()
        {
            mesh.UnperformSetZero(edge, newState, edgeSetChanges);
        }

#endregion

#region IEquatable<IAction> Members

        public bool Equals(IAction other)
        {
            SetAction otherAction = other as SetAction;
            if (otherAction == null)
                return false;
            else
                return otherAction.newState == this.newState && otherAction.mesh == this.mesh && otherAction.edge == this.edge;
        }

#endregion

        public override int GetHashCode()
        {
            return HashCode.Combine(mesh.GetHashCode(), edge.GetHashCode(), newState.GetHashCode());
        }

        public override bool Equals(object obj)
        {
            if (obj is IAction)
                return Equals((IAction)obj);
            return false;
        }
    }

    public class UnsetAction : IAction
    {

        public UnsetAction(Mesh mesh, int edge)
        {
            this.mesh = mesh;
            this.edge = edge;
        }

        private Mesh mesh;

        public int EdgeIndex
        {
            get
            {
                return edge;
            }
        }
        private int edge;
        private EdgeState oldState;
        private List<int[]> edgeSetChanges;

        public bool Successful
        {
            get
            {
                return successful;
            }
        }
        private bool successful;

#region IAction Members

        public string Name
        {
            get { return "Clear Edge " + edge.ToString(); }
        }

        public bool Perform()
        {
            successful = true;
            edgeSetChanges = new List<int[]>();
            if (!mesh.PerformUnsetZero(edge, out oldState, edgeSetChanges))
            {
                successful = false;
            }
            return true;
        }

        public void Unperform()
        {
            mesh.UnperformUnsetZero(edge, oldState, edgeSetChanges);
        }

#endregion

#region IEquatable<IAction> Members

        public bool Equals(IAction other)
        {
            UnsetAction otherAction = other as UnsetAction;
            if (otherAction == null)
                return false;
            else
                return otherAction.mesh == this.mesh && otherAction.edge == this.edge;
        }

#endregion

        public override int GetHashCode()
        {
            return HashCode.Combine(mesh.GetHashCode(), edge.GetHashCode());
        }

        public override bool Equals(object obj)
        {
            if (obj is IAction)
                return Equals((IAction)obj);
            return false;
        }
    }

    public class ColorJoinAction : IAction
    {

        public ColorJoinAction(Mesh mesh, int edge1, int edge2, bool same)
        {
            this.mesh = mesh;
            this.edge1 = edge1;
            this.edge2 = edge2;
            this.same = same;
        }

        private Mesh mesh;
        public int Edge1
        {
            get
            {
                return edge1;
            }
        }
        private int edge1;
        public int Edge2
        {
            get
            {
                return edge2;
            }
        }
        private int edge2;
        public bool Same
        {
            get
            {
                return same;
            }
        }
        private bool same;

        public List<int> GetAffectedEdges()
        {
            List<int> res = new List<int>();
            foreach (int[] change in colorSetChanges)
            {
                // Ignore changes to the number of sets.
                if (change.Length > 1)
                    res.Add(change[0]);
            }
            return res;
        }

        private List<int[]> colorSetChanges;

        public bool Successful
        {
            get
            {
                return successful;
            }
        }
        private bool successful;

        public bool WasteOfTime
        {
            get
            {
                return wasteOfTime;
            }
        }
        private bool wasteOfTime;

#region IAction Members

        public string Name
        {
            get { return "Color Edges " + edge1.ToString() + " and " + edge2.ToString() + (same ? " the same" : " opposite"); }
        }

        public bool Perform()
        {
            successful = true;
            colorSetChanges = new List<int[]>();
            if (!mesh.JoinColor(edge1, edge2, !same, colorSetChanges, ref wasteOfTime))
            {
                return false;
            }
            return true;
        }

        public void Unperform()
        {
            mesh.UnjoinColor(colorSetChanges);
        }

#endregion

#region IEquatable<IAction> Members

        public bool Equals(IAction other)
        {
            ColorJoinAction otherAction = other as ColorJoinAction;
            if (otherAction == null)
                return false;
            else
                return otherAction.mesh == this.mesh && otherAction.edge1 == this.edge1 && otherAction.edge2 == this.edge2 && otherAction.same == this.same;
        }

#endregion

        public override int GetHashCode()
        {
            return HashCode.Combine(mesh.GetHashCode(), edge1.GetHashCode(), edge2.GetHashCode(), same.GetHashCode());
        }

        public override bool Equals(object obj)
        {
            if (obj is IAction)
                return Equals((IAction)obj);
            return false;
        }
    }

    public class CellColorJoinAction : IAction
    {

        public CellColorJoinAction(Mesh mesh, int cell1, int cell2, bool same)
        {
            this.mesh = mesh;
            this.cell1 = cell1;
            this.cell2 = cell2;
            this.same = same;
#if DEBUG
            //st = new StackTrace(true);
#endif
        }
#if DEBUG
        //public StackTrace st;
#endif
        private Mesh mesh;
        public int Cell1
        {
            get
            {
                return cell1;
            }
        }
        private int cell1;
        public int Cell2
        {
            get
            {
                return cell2;
            }
        }
        private int cell2;
        public bool Same
        {
            get
            {
                return same;
            }
        }
        private bool same;

        public List<int> GetAffectedCells()
        {
            List<int> res = new List<int>();
            foreach (int[] change in colorSetChanges)
            {
                // Ignore changes to the number of sets.
                if (change.Length > 1)
                    res.Add(change[0]);
            }
            return res;
        }

        private List<int[]> colorSetChanges;

        public bool Successful
        {
            get
            {
                return successful;
            }
        }
        private bool successful;

        public bool WasteOfTime
        {
            get
            {
                return wasteOfTime;
            }
        }
        private bool wasteOfTime;

#region IAction Members

        public string Name
        {
            get { return "Color Cells " + cell1.ToString() + " and " + cell2.ToString() + (same ? " the same" : " opposite"); }
        }

        public bool Perform()
        {
            successful = true;
            colorSetChanges = new List<int[]>();
            if (!mesh.JoinCellColor(cell1, cell2, !same, colorSetChanges, ref wasteOfTime))
            {
                return false;
            }
            return true;
        }

        public void Unperform()
        {
            mesh.UnjoinCellColor(colorSetChanges);
        }

#endregion

#region IEquatable<IAction> Members

        public bool Equals(IAction other)
        {
            CellColorJoinAction otherAction = other as CellColorJoinAction;
            if (otherAction == null)
                return false;
            else
                return otherAction.mesh == this.mesh && otherAction.cell1 == this.cell1 && otherAction.cell2 == this.cell2 && otherAction.same == this.same;
        }

#endregion

        public override int GetHashCode()
        {
            return HashCode.Combine(mesh.GetHashCode(), cell1.GetHashCode(), cell1.GetHashCode(), same.GetHashCode());
        }

        public override bool Equals(object obj)
        {
            if (obj is IAction)
                return Equals((IAction)obj);
            return false;
        }
    }

    public class CellColorClearAction : IAction
    {

        public CellColorClearAction(Mesh mesh, int cell1)
        {
            this.mesh = mesh;
            this.cell1 = cell1;
        }

        private Mesh mesh;
        public int Cell1
        {
            get
            {
                return cell1;
            }
        }
        private int cell1;

        private int oldColor;
        private int oldSetPos;


        public bool Successful
        {
            get
            {
                return successful;
            }
        }
        private bool successful;

#region IAction Members

        public string Name
        {
            get { return "Clear Color of Cell " + cell1.ToString(); }
        }

        public bool Perform()
        {
            successful = true;
            oldColor = mesh.Cells[cell1].Color;
            if (!mesh.ClearCellColor(cell1, out oldSetPos))
            {
                return false;
            }
            return true;
        }

        public void Unperform()
        {
            mesh.ResetCellColor(cell1, oldColor, oldSetPos);
        }

#endregion

#region IEquatable<IAction> Members

        public bool Equals(IAction other)
        {
            CellColorClearAction otherAction = other as CellColorClearAction;
            if (otherAction == null)
                return false;
            else
                return otherAction.mesh == this.mesh && otherAction.cell1 == this.cell1;
        }

#endregion

        public override int GetHashCode()
        {
            return HashCode.Combine(mesh.GetHashCode(), cell1.GetHashCode());
        }

        public override bool Equals(object obj)
        {
            if (obj is IAction)
                return Equals((IAction)obj);
            return false;
        }
    }

    public class EdgeRestrictionAction : IAction
    {

        public EdgeRestrictionAction(Mesh mesh, int edge1, int edge2, EdgePairRestriction state)
        {
            this.mesh = mesh;
            this.edge1 = edge1;
            this.edge2 = edge2;
            this.state = state;
#if DEBUG
            //st = new StackTrace(true);
#endif
        }
#if DEBUG
        //public StackTrace st;
#endif
        private Mesh mesh;
        public int Edge1
        {
            get
            {
                return edge1;
            }
        }
        private int edge1;
        public int Edge2
        {
            get
            {
                return edge2;
            }
        }
        private int edge2;
        public EdgePairRestriction State
        {
            get
            {
                return state;
            }
        }
        private EdgePairRestriction state;

        public List<int> GetAffectedEdges()
        {
            List<int> res = new List<int>();
            if (!wasteOfTime)
            {
                res.Add(edge1);
                res.Add(edge2);
            }
            return res;
        }

        public bool Successful
        {
            get
            {
                return successful;
            }
        }
        private bool successful;

        public bool WasteOfTime
        {
            get
            {
                return wasteOfTime;
            }
        }
        private bool wasteOfTime;

#region IAction Members

        public string Name
        {
            get { return "Restrict Edges " + edge1.ToString() + " and " + edge2.ToString() + " to be " + state.ToString(); }
        }

        public bool Perform()
        {
            successful = true;
            if (!mesh.RestrictEdges(edge1, edge2, state, ref wasteOfTime))
            {
                return false;
            }
            return true;
        }

        public void Unperform()
        {
            if (!wasteOfTime)
                mesh.UnRestrictEdges(edge1, edge2);
        }

#endregion

#region IEquatable<IAction> Members

        public bool Equals(IAction other)
        {
            EdgeRestrictionAction otherAction = other as EdgeRestrictionAction;
            if (otherAction == null)
                return false;
            else
                return otherAction.mesh == this.mesh && otherAction.edge1 == this.edge1 && otherAction.edge2 == this.edge2 && otherAction.state == this.state;
        }

#endregion

        public override int GetHashCode()
        {
            return HashCode.Combine(mesh.GetHashCode(), edge1.GetHashCode(), edge2.GetHashCode(), ((int)state).GetHashCode());
        }

        public override bool Equals(object obj)
        {
            if (obj is IAction)
                return Equals((IAction)obj);
            return false;
        }
    }

    public static class HashCode
    {
        public static int Combine(int a, int b)
        {
            return ((a << 5) + a) ^ b;
        }
        public static int Combine(int a, int b, int c)
        {
            int middle= ((a << 5) + a) ^ b;
            return ((middle << 5) + middle) ^ c;
        }
        public static int Combine(int a, int b, int c, int d)
        {
            int middle = ((a << 5) + a) ^ b;
            middle = ((middle << 5) + middle) ^ c;
            return ((middle << 5) + middle) ^ d;
        }
    }

    public class ActionSorter : IComparer<IAction>
    {
#region IComparer<IAction> Members

        public int Compare(IAction x, IAction y)
        {
            if (x is SetAction)
            {
                if (y is SetAction)
                {
                    return ((SetAction)x).EdgeIndex.CompareTo(((SetAction)y).EdgeIndex);
                }
                else
                {
                    return 1;
                }
            }
            else if (x is ColorJoinAction)
            {
                ColorJoinAction xc = (ColorJoinAction)x;
                if (y is SetAction)
                {
                    return -1;
                }
                else if (y is ColorJoinAction)
                {
                    ColorJoinAction yc = (ColorJoinAction)y;
                    int res = xc.Edge1.CompareTo(yc.Edge1);
                    if (res == 0)
                        res = xc.Edge2.CompareTo(yc.Edge2);
                    return res;
                }
                else
                {
                    return 1;
                }
            }
            else if (x is CellColorJoinAction)
            {
                CellColorJoinAction xc = (CellColorJoinAction)x;
                if (y is CellColorJoinAction)
                {
                    CellColorJoinAction yc = (CellColorJoinAction)y;
                    int res = xc.Cell1.CompareTo(yc.Cell1);
                    if (res == 0)
                        res = xc.Cell2.CompareTo(yc.Cell2);
                    return res;

                }
                else
                {
                    return -1;
                }
            }
            throw new NotSupportedException("ActionSorter doesn't support the provided action type.");
        }

#endregion
    }

    public class Chain
    {
        public ChainNode Start;
        public ChainNode End;
    }

    public class ChainNode
    {
        public int Intersection;

        public ChainNode Link1;
        public ChainNode Link2;

        public ChainNode GetNext(ChainNode prior)
        {
            if (prior == Link1)
                return Link2;
            return Link1;
        }

        public void Join(ChainNode other)
        {
            if (Link1 == null)
            {
                Link1 = other;
                Link1.JoinInternal(this);
            }
            else if (Link2 == null)
            {
                Link2 = other;
                Link2.JoinInternal(this);
            }
            else
                throw new InvalidOperationException("Can't join chain, already doubly joined.");
        }
        private void JoinInternal(ChainNode other)
        {
            if (Link1 == null)
            {
                Link1 = other;
            }
            else if (Link2 == null)
            {
                Link2 = other;
            }
            else
                throw new InvalidOperationException("Can't join chain, already doubly joined.");
        }

        public void Remove(ChainNode other)
        {
            if (Link1 == other)
            {
                Link1 = null;
                other.RemoveInternal(this);
            }
            else if (Link2 == other)
            {
                Link2 = null;
                other.RemoveInternal(this);
            }
            else
                throw new InvalidOperationException("Can't break from chain, not joined.");
        }
        private void RemoveInternal(ChainNode other)
        {
            if (Link1 == other)
            {
                Link1 = null;
            }
            else if (Link2 == other)
            {
                Link2 = null;
            }
            else
                throw new InvalidOperationException("Can't break from chain, not joined.");
        }

        public int[] EdgesToLink1;
        public int[] EdgesToLink2;
    }

    /// <summary>
    /// Tracks which set each item belongs to after a set of union operations to merge disjoint sets.
    /// This is very efficient, it can be consider that each operation takes close to O(1) time on average.
    /// </summary>
    public class DisjointTracker
    {

        /// <summary>
        /// Constructor.
        /// </summary>
        public DisjointTracker(int size)
        {
            this.tracker = new List<int>(size);
            this.ranker = new List<int>(size);
        }

        private List<int> tracker;
        private List<int> ranker;


        /// <summary>
        /// Adds a new item to be tracked.  Initially it is in a set of its own.
        /// </summary>
        /// <param name="value">
        /// Value to be added to the tracking.
        /// </param>
        public void Add(int value)
        {
            ResizeIfNeeded(value);
            tracker[value] = value;
            ranker[value] = 0;
        }

        private void ResizeIfNeeded(int value)
        {
            while (value >= tracker.Count)
            {
                tracker.Add(-1);
                ranker.Add(-1);
            }
        }

        /// <summary>
        /// Unions the two sets which contain the specified items.  Does nothing if the items are already in the same set.
        /// </summary>
        /// <param name="first">
        /// First item to find set to union.
        /// </param>
        /// <param name="second">
        /// Second item to find set to union.
        /// </param>
        public void Union(int first, int second)
        {
            Link(GetRepresentative(first), GetRepresentative(second));
        }

        /// <summary>
        /// Gets the primary member of the set which a given value belongs to.
        /// </summary>
        /// <param name="value">
        /// Value to lookup.
        /// </param>
        /// <returns>
        /// The primary member of the set that the value is currently a member of.
        /// </returns>
        public int GetRepresentative(int value)
        {
            int parent = tracker[value];
            if (parent == value)
                return parent;
            int realParent = GetRepresentative(parent);
            tracker[value] = realParent;
            return realParent;
        }

        private void Link(int first, int second)
        {
            if (first == second)
                return;
            int firstRank = ranker[first];
            int secondRank = ranker[second];
            if (firstRank > secondRank)
            {
                tracker[second] = first;
            }
            else
            {
                if (firstRank == secondRank)
                    ranker[second] = secondRank + 1;
                tracker[first] = second;
            }
        }
    }

    // Dirty evil list replacement.  Minimal functionality and convienience, maximal speed?
    public class QuickList
    {
        public QuickList()
            : this(4)
        {
        }
        public QuickList(int length)
        {
            Buffer = new int[length];
        }
        public int Count;
        public void Add(int value)
        {
            if (Count == Buffer.Length)
            {
                int[] newBuffer = new int[Buffer.Length * 2];
                Array.Copy(Buffer, newBuffer, Count);
                Buffer = newBuffer;
            }
            Buffer[Count] = value;
            Count++;
        }
        public void Clear()
        {
            Count = 0;
        }

        public int[] Buffer;
    }

}
