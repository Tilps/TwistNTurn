using System;
using System.Collections.Generic;

namespace TwistNTurn
{
    public partial class Mesh
    {
        private void ConstructGrid(int width, int height, MeshType type)
        {
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
                            for (var index1 = 0; index1 < newEdge.Cells.Count; index1++)
                            {
                                int cellIndex = newEdge.Cells[index1];
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
                            for (var index1 = 0; index1 < newEdge.Cells.Count; index1++)
                            {
                                int cellIndex = newEdge.Cells[index1];
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
                        for (int j = 0; j < cellCorners.Count; j++)
                        {
                            int cellCorner = cellCorners[j];
                            cell.Intersections.Add(cellCorner);
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
            for (var index = 0; index < inters.Edges.Count; index++)
            {
                int eIndex = inters.Edges[index];
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
                    for (var index1 = 0; index1 < nextInter.Edges.Count; index1++)
                    {
                        int edgeIndex = nextInter.Edges[index1];
                        int possibleNext = GetOtherInters(edgeIndex, last);
                        if (possibleNext == realLast)
                            continue;
                        Intersection possibleNextInter = intersections[possibleNext];
                        double possibleNextAngle =
                            Math.Atan2(possibleNextInter.Y - lastInter.Y, possibleNextInter.X - lastInter.X);
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
                        double angle = Math.Acos(dot / Math.Sqrt(dx1 * dx1 + dy1 * dy1) /
                                                 Math.Sqrt(dx2 * dx2 + dy2 * dy2));
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

    }
}