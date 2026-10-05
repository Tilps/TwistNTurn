using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

namespace TwistNTurn
{
    public partial class Mesh
    {
        public bool LoadFromText(ReadOnlySpan<char> text)
        {
            List<string> lines = new List<string>();
            foreach (var line in text.EnumerateLines())
            {
                lines.Add(line.ToString());
            }
            return LoadFromText(lines.ToArray());
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
            if (edgePairRestrictions != null)
            {
                for (int i = 0; i < edges.Count; i++)
                    for (int j = i + 1; j < edges.Count; j++)
                        if (edgePairRestrictions[i, j] != EdgePairRestriction.None)
                            counter++;
            }
            writer.WriteLine(counter);
            if (edgePairRestrictions != null)
            {
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
        }

        public string SaveToString()
        {
            using var sw = new StringWriter();
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
            int width;
            int height;
            if (!int.TryParse(dimensions[0], out width) || width < 2)
                throw new ArgumentException("puzzle string invalid.");
            if (!int.TryParse(dimensions[1], out height) || height < 2)
                throw new ArgumentException("puzzle string invalid.");
            Mesh result = new Mesh(width - 1, height - 1, MeshType.Square);
            if (bits.Length > 1)
            {
                string clueStr = bits[1];
                if (clueStr.Contains('.') || clueStr.Contains('2'))
                {
                    // Sequential intersection clue string (e.g. '.' for unknown, '1' for twist, '2' for turn)
                    for (int i = 0; i < clueStr.Length && i < result.Intersections.Count; i++)
                    {
                        char c = clueStr[i];
                        if (c == '1' || c == 'w' || c == 'W')
                            result.Intersections[i].Type = IntersType.Twist;
                        else if (c == '2' || c == 'b' || c == 'B')
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
                        else if (c >= 'a' && c <= 'z')
                        {
                            counter += (int)(c - 'a' + 1);
                        }
                        else
                        {
                            counter++;
                        }
                    }
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
                    result.Append(count);
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
                    currentCount = currentCount * 10 + (c - '0');
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