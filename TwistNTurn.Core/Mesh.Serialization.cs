using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

namespace TwistNTurn
{
    public partial class Mesh
    {
        private const int MaxSerializedIntersections = 50_000;
        private const int MaxSerializedEdges = 50_000;

        public bool LoadFromText(ReadOnlySpan<char> text)
        {
            List<string> lines = [];
            foreach (var line in text.EnumerateLines())
            {
                lines.Add(line.ToString());
            }
            return LoadFromText(lines.ToArray());
        }

        private static bool TryReadLine(string[] lines, ref int curLine, out ReadOnlySpan<char> line)
        {
            if (curLine < 0 || curLine >= lines.Length)
            {
                line = default;
                return false;
            }
            line = lines[curLine++].AsSpan().Trim();
            return true;
        }

        private static bool TryReadInt(string[] lines, ref int curLine, out int value)
        {
            value = 0;
            return TryReadLine(lines, ref curLine, out var span)
                && int.TryParse(span, NumberStyles.Integer, CultureInfo.InvariantCulture, out value);
        }

        public bool LoadFromText(string[] lines)
        {
            if (lines == null || lines.Length < 4)
                return false;

            ReadOnlySpan<char> header = lines[0].AsSpan().Trim();
            if (!Enum.TryParse<MeshType>(header, ignoreCase: true, out MeshType parsedType) || !Enum.IsDefined(parsedType))
                return false;

            int curLine = 1;
            if (!TryReadLine(lines, ref curLine, out var section) || !section.SequenceEqual("Intersections"))
                return false;
            if (!TryReadInt(lines, ref curLine, out int intersectionCount) || intersectionCount < 0 || intersectionCount > MaxSerializedIntersections)
                return false;

            var parsedIntersections = new List<Intersection>(intersectionCount);
            for (int i = 0; i < intersectionCount; i++)
            {
                if (!TryReadLine(lines, ref curLine, out var span))
                    return false;
                string[] splits = span.ToString().Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (splits.Length < 3)
                    return false;
                if (!float.TryParse(splits[0], NumberStyles.Float, CultureInfo.InvariantCulture, out float x) || !float.IsFinite(x))
                    return false;
                if (!float.TryParse(splits[1], NumberStyles.Float, CultureInfo.InvariantCulture, out float y) || !float.IsFinite(y))
                    return false;
                if (!Enum.TryParse<IntersType>(splits[2], ignoreCase: true, out IntersType interType) || !Enum.IsDefined(interType))
                    return false;

                parsedIntersections.Add(new Intersection { X = x, Y = y, Type = interType });
            }

            if (!TryReadLine(lines, ref curLine, out section) || !section.SequenceEqual("Edges"))
                return false;
            if (!TryReadInt(lines, ref curLine, out int edgeCount) || edgeCount < 0 || edgeCount > MaxSerializedEdges)
                return false;

            var parsedEdges = new List<(int V1, int V2, EdgeState State)>(edgeCount);
            var seenEdgePairs = new HashSet<(int, int)>(edgeCount);
            for (int i = 0; i < edgeCount; i++)
            {
                if (!TryReadLine(lines, ref curLine, out var span))
                    return false;
                string[] splits = span.ToString().Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (splits.Length < 3)
                    return false;
                if (!int.TryParse(splits[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out int v1)
                    || !int.TryParse(splits[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int v2))
                    return false;
                if (v1 < 0 || v1 >= intersectionCount || v2 < 0 || v2 >= intersectionCount || v1 == v2)
                    return false;
                if (!Enum.TryParse<EdgeState>(splits[2], ignoreCase: true, out EdgeState state) || !Enum.IsDefined(state))
                    return false;

                var key = v1 < v2 ? (v1, v2) : (v2, v1);
                if (!seenEdgePairs.Add(key))
                    return false;

                parsedEdges.Add((v1, v2, state));
            }

            // Parse optional sections into actions after mesh topology is built; validate structure first.
            int optionalStartLine = curLine;
            var edgeColorJoins = new List<(int First, int Edge, bool Same)>();
            var cellColorJoins = new List<(int Cell1, int Cell2, bool Same)>();
            var edgeRestrictions = new List<(int Edge1, int Edge2, EdgePairRestriction Restriction)>();

            if (curLine < lines.Length && lines[curLine].Trim() == "EdgeColorSets")
            {
                curLine++;
                if (!TryReadInt(lines, ref curLine, out int edgeColorSetsCount) || edgeColorSetsCount < 0 || edgeColorSetsCount > MaxSerializedEdges)
                    return false;
                for (int i = 0; i < edgeColorSetsCount; i++)
                {
                    if (!TryReadLine(lines, ref curLine, out _)) return false; // "EdgeColorSet"
                    if (!TryReadInt(lines, ref curLine, out int edgeColorSetCount) || edgeColorSetCount < 0 || edgeColorSetCount > edgeCount)
                        return false;
                    int first = -1;
                    int firstColor = 0;
                    for (int j = 0; j < edgeColorSetCount; j++)
                    {
                        if (!TryReadLine(lines, ref curLine, out var lineSpan)) return false;
                        string[] bits = lineSpan.ToString().Split(' ', StringSplitOptions.RemoveEmptyEntries);
                        if (bits.Length < 2
                            || !int.TryParse(bits[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out int edgeIdx)
                            || !int.TryParse(bits[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int color)
                            || edgeIdx < 0 || edgeIdx >= edgeCount)
                            return false;

                        if (first == -1)
                        {
                            first = edgeIdx;
                            firstColor = color;
                        }
                        else
                        {
                            edgeColorJoins.Add((first, edgeIdx, color == firstColor));
                        }
                    }
                }
            }

            if (curLine < lines.Length && lines[curLine].Trim() == "CellColorSets")
            {
                curLine++;
                if (!TryReadInt(lines, ref curLine, out int cellColorSetsCount) || cellColorSetsCount < 0 || cellColorSetsCount > MaxSerializedEdges)
                    return false;
                for (int i = 0; i < cellColorSetsCount; i++)
                {
                    if (!TryReadLine(lines, ref curLine, out _)) return false; // "CellColorSet"
                    if (!TryReadInt(lines, ref curLine, out int cellColorSetCount) || cellColorSetCount < 0 || cellColorSetCount > MaxSerializedEdges)
                        return false;
                    int first = -1;
                    int firstColor = 0;
                    for (int j = 0; j < cellColorSetCount; j++)
                    {
                        if (!TryReadLine(lines, ref curLine, out var lineSpan)) return false;
                        string[] bits = lineSpan.ToString().Split(' ', StringSplitOptions.RemoveEmptyEntries);
                        if (bits.Length < 2
                            || !int.TryParse(bits[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out int cellIdx)
                            || !int.TryParse(bits[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int color)
                            || cellIdx < 0)
                            return false;

                        if (Math.Abs(color) != 1)
                        {
                            if (first == -1)
                            {
                                first = cellIdx;
                                firstColor = color;
                            }
                            else
                            {
                                cellColorJoins.Add((first, cellIdx, color == firstColor));
                            }
                        }
                        else
                        {
                            cellColorJoins.Add((cellIdx, -1, color == 1));
                        }
                    }
                }
            }

            if (curLine < lines.Length && lines[curLine].Trim() == "EdgePairRestrictions")
            {
                curLine++;
                if (!TryReadInt(lines, ref curLine, out int edgePairRestrictionCount) || edgePairRestrictionCount < 0 || edgePairRestrictionCount > MaxSerializedEdges * 4)
                    return false;
                for (int j = 0; j < edgePairRestrictionCount; j++)
                {
                    if (!TryReadLine(lines, ref curLine, out var lineSpan)) return false;
                    string[] bits = lineSpan.ToString().Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    if (bits.Length < 3
                        || !int.TryParse(bits[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out int edge1)
                        || !int.TryParse(bits[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int edge2)
                        || edge1 < 0 || edge1 >= edgeCount || edge2 < 0 || edge2 >= edgeCount
                        || !Enum.TryParse<EdgePairRestriction>(bits[2], ignoreCase: true, out EdgePairRestriction restriction)
                        || !Enum.IsDefined(restriction))
                        return false;

                    edgeRestrictions.Add((edge1, edge2, restriction));
                }
            }

            _ = optionalStartLine;

            // Commit validated mesh atomically
            meshType = parsedType;
            Intersections.Clear();
            Edges.Clear();
            Cells.Clear();
            edgeDistancesCache = null;
            boringEdges = null;
            patternLookup.Clear();

            Intersections.AddRange(parsedIntersections);
            List<IAction> settings = [];
            for (int i = 0; i < parsedEdges.Count; i++)
            {
                var (v1, v2, state) = parsedEdges[i];
                AddEdge(v1, v2);
                if (state != EdgeState.Empty)
                {
                    settings.Add(new SetAction(this, i, state));
                }
            }

            CreateOpposites();
            CreateCells();
            Clear();

            foreach (var (first, edgeIdx, same) in edgeColorJoins)
                settings.Add(new ColorJoinAction(this, first, edgeIdx, same));
            foreach (var (c1, c2, same) in cellColorJoins)
            {
                if (c1 < Cells.Count && (c2 == -1 || c2 < Cells.Count))
                    settings.Add(new CellColorJoinAction(this, c1, c2, same));
            }
            foreach (var (e1, e2, restriction) in edgeRestrictions)
                settings.Add(new EdgeRestrictionAction(this, e1, e2, restriction));

            PerformListRegardless(settings);
            for (int i = 0; i < intersections.Count; i++)
            {
                if (intersections[i].Edges.Count != intersections[i].Opposite.Count)
                    return false;
            }
            return true;
        }

        public void Save(TextWriter writer)
        {
            writer.WriteLine(MeshType.ToString());
            writer.WriteLine("Intersections");
            writer.WriteLine(Intersections.Count.ToString(CultureInfo.InvariantCulture));
            foreach (Intersection inters in Intersections)
            {
                writer.Write(inters.X.ToString(CultureInfo.InvariantCulture));
                writer.Write(' ');
                writer.Write(inters.Y.ToString(CultureInfo.InvariantCulture));
                writer.Write(' ');
                writer.Write(inters.Type.ToString());
                writer.WriteLine();
            }
            writer.WriteLine("Edges");
            writer.WriteLine(Edges.Count.ToString(CultureInfo.InvariantCulture));
            foreach (Edge edge in Edges)
            {
                writer.Write(edge.Intersections[0].ToString(CultureInfo.InvariantCulture));
                writer.Write(' ');
                writer.Write(edge.Intersections[1].ToString(CultureInfo.InvariantCulture));
                writer.Write(' ');
                writer.WriteLine(edge.State.ToString());
            }
            writer.WriteLine("EdgeColorSets");
            writer.WriteLine(colorSets.Count.ToString(CultureInfo.InvariantCulture));
            foreach (List<int> colorSet in colorSets)
            {
                writer.WriteLine("EdgeColorSet");
                writer.WriteLine(colorSet.Count.ToString(CultureInfo.InvariantCulture));
                foreach (int edge in colorSet)
                {
                    writer.Write(edge.ToString(CultureInfo.InvariantCulture));
                    writer.Write(' ');
                    writer.WriteLine(edges[edge].Color.ToString(CultureInfo.InvariantCulture));
                }
            }
            writer.WriteLine("CellColorSets");
            writer.WriteLine(cellColorSets.Count.ToString(CultureInfo.InvariantCulture));
            foreach (List<int> colorSet in cellColorSets)
            {
                writer.WriteLine("CellColorSet");
                writer.WriteLine(colorSet.Count.ToString(CultureInfo.InvariantCulture));
                foreach (int cell in colorSet)
                {
                    writer.Write(cell.ToString(CultureInfo.InvariantCulture));
                    writer.Write(' ');
                    writer.WriteLine(cells[cell].Color.ToString(CultureInfo.InvariantCulture));
                }
            }
            writer.WriteLine("EdgePairRestrictions");
            int counter = 0;
            if (edgePairRestrictions != null)
            {
                for (int i = 0; i < edges.Count; i++)
                    for (int j = i + 1; j < edges.Count; j++)
                        if (edgePairRestrictions[i, j] != EdgePairRestriction.None)
                            counter++;
            }
            writer.WriteLine(counter.ToString(CultureInfo.InvariantCulture));
            if (edgePairRestrictions != null)
            {
                for (int i = 0; i < edges.Count; i++)
                    for (int j = i + 1; j < edges.Count; j++)
                        if (edgePairRestrictions[i, j] != EdgePairRestriction.None)
                        {
                            writer.Write(i.ToString(CultureInfo.InvariantCulture));
                            writer.Write(' ');
                            writer.Write(j.ToString(CultureInfo.InvariantCulture));
                            writer.Write(' ');
                            writer.WriteLine(edgePairRestrictions[i, j].ToString());
                        }
            }
        }

        public string SaveToString()
        {
            using var sw = new StringWriter(CultureInfo.InvariantCulture);
            Save(sw);
            return sw.ToString();
        }

        internal static Mesh FromPuzzleString(string puzzle)
        {
            string[] bits = puzzle.Split(':');
            if (bits.Length != 2)
                throw new ArgumentException("puzzle string invalid.");
            string[] dimensions = bits[0].Split('x');
            if (dimensions.Length != 2)
                throw new ArgumentException("puzzle string invalid.");
            if (!int.TryParse(dimensions[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out int width)
                || width < 2 || width - 1 > PuzzleHelper.MaxGridDimension)
                throw new ArgumentException("puzzle string invalid.");
            if (!int.TryParse(dimensions[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int height)
                || height < 2 || height - 1 > PuzzleHelper.MaxGridDimension)
                throw new ArgumentException("puzzle string invalid.");

            Mesh result = new Mesh(width - 1, height - 1, MeshType.Square);
            string clueStr = bits[1];
            if (clueStr.Contains('.') || clueStr.Contains('2'))
            {
                for (int i = 0; i < clueStr.Length && i < result.Intersections.Count; i++)
                {
                    char c = clueStr[i];
                    if (c is '1' or 'w' or 'W')
                        result.Intersections[i].Type = IntersType.Twist;
                    else if (c is '2' or 'b' or 'B')
                        result.Intersections[i].Type = IntersType.Turn;
                    else
                        result.Intersections[i].Type = IntersType.Unknown;
                }
            }
            else
            {
                int counter = 0;
                foreach (char c in clueStr)
                {
                    if (counter >= width * height) break;
                    int x = counter % width;
                    int y = counter / width;
                    int pos = y + x * height;
                    if (c == '0')
                    {
                        if (pos < result.Intersections.Count)
                            result.Intersections[pos].Type = IntersType.Twist;
                        counter++;
                    }
                    else if (c == '1')
                    {
                        if (pos < result.Intersections.Count)
                            result.Intersections[pos].Type = IntersType.Turn;
                        counter++;
                    }
                    else if (c is >= 'a' and <= 'z')
                    {
                        counter += (c - 'a' + 1);
                    }
                    else
                    {
                        counter++;
                    }
                }
            }
            return result;
        }

        internal string ToPuzzleString()
        {
            if (!PuzzleFormatHelper.TryGetSquareDimensions(this, out int cellW, out int cellH))
                return string.Empty;

            int width = cellW + 1;
            int height = cellH + 1;
            StringBuilder result = new StringBuilder();
            result.AppendFormat(CultureInfo.InvariantCulture, "{0}x{1}:", width, height);
            int counter = 0;
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int pos = y + x * height;
                    switch (intersections[pos].Type)
                    {
                        case IntersType.Twist:
                            if (counter > 0)
                                result.Append((char)('a' + counter - 1));
                            result.Append('0');
                            counter = 0;
                            break;
                        case IntersType.Turn:
                            if (counter > 0)
                                result.Append((char)('a' + counter - 1));
                            result.Append('1');
                            counter = 0;
                            break;
                        case IntersType.Unknown:
                            counter++;
                            if (counter > 23)
                            {
                                result.Append((char)('a' + counter - 1));
                                counter = 0;
                            }
                            break;
                    }
                }
            }
            return result.ToString();
        }

        internal string ToEdgeStateString()
        {
            if (edges.Count == 0) return string.Empty;
            StringBuilder result = new StringBuilder();
            int i = 0;
            while (i < edges.Count)
            {
                EdgeState state = edges[i].State;
                char code = state switch
                {
                    EdgeState.Filled => 'F',
                    EdgeState.Excluded => 'X',
                    _ => '.'
                };
                int count = 1;
                while (i + count < edges.Count && edges[i + count].State == state)
                {
                    count++;
                }
                if (count > 1)
                {
                    result.Append(count.ToString(CultureInfo.InvariantCulture));
                }
                result.Append(code);
                i += count;
            }
            return result.ToString();
        }

        internal void FromEdgeStateString(string input)
        {
            Clear();
            if (string.IsNullOrEmpty(input))
            {
                RecalculateCounts();
                return;
            }

            int edgeIdx = 0;
            int currentCount = 0;
            for (int i = 0; i < input.Length; i++)
            {
                char c = input[i];
                if (char.IsDigit(c))
                {
                    if (currentCount <= edges.Count)
                        currentCount = Math.Min(edges.Count + 1, currentCount * 10 + (c - '0'));
                }
                else
                {
                    int count = currentCount == 0 ? 1 : currentCount;
                    currentCount = 0;

                    EdgeState state = c switch
                    {
                        'F' or '1' => EdgeState.Filled,
                        'X' or 'x' or '2' => EdgeState.Excluded,
                        _ => EdgeState.Empty
                    };

                    for (int k = 0; k < count && edgeIdx < edges.Count; k++)
                    {
                        edges[edgeIdx++].State = state;
                    }
                }
            }
            RecalculateCounts();
        }
    }
}
