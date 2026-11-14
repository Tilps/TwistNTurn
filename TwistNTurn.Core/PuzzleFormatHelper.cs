using System;
using System.Text;

namespace TwistNTurn
{
    /// <summary>
    /// Converts puzzles between native .twist text, compact share codes, pearl descriptors, and Masyu ASCII grids.
    /// </summary>
    public static class PuzzleFormatHelper
    {
        private static readonly char[] LineSeparators = ['\r', '\n'];

        public static bool CanExportSquareFormats(Mesh? mesh)
            => mesh != null && TryGetSquareDimensions(mesh, out _, out _);

        public static bool TryGetSquareDimensions(Mesh mesh, out int width, out int height)
        {
            width = 0;
            height = 0;
            if (mesh.MeshType != MeshType.Square && mesh.MeshType != MeshType.SquareSymmetrical)
                return false;

            if (mesh.Cells.Count == 0 || mesh.Intersections.Count < 4)
                return false;

            int v = 0;
            while (v + 1 < mesh.Intersections.Count && TryFindEdgeJoining(mesh, v, v + 1, out _))
            {
                v++;
            }

            if (v <= 0 || mesh.Intersections.Count % (v + 1) != 0)
                return false;

            height = v;
            width = (mesh.Intersections.Count / (height + 1)) - 1;
            return width > 0 && mesh.Cells.Count == width * height;
        }

        private static bool TryFindEdgeJoining(Mesh mesh, int v1, int v2, out int edgeIndex)
        {
            edgeIndex = -1;
            if (v1 < 0 || v2 < 0 || v1 >= mesh.Intersections.Count || v2 >= mesh.Intersections.Count)
                return false;

            var e1 = mesh.Intersections[v1].Edges;
            var e2 = mesh.Intersections[v2].Edges;
            for (int i = 0; i < e1.Count; i++)
            {
                int candidate = e1[i];
                for (int j = 0; j < e2.Count; j++)
                {
                    if (e2[j] == candidate)
                    {
                        edgeIndex = candidate;
                        return true;
                    }
                }
            }
            return false;
        }

        public static string ExportTwistText(Mesh mesh) => mesh.SaveToString();

        public static string? ExportPearlDescriptor(Mesh mesh)
        {
            if (!TryGetSquareDimensions(mesh, out _, out _))
                return null;
            return mesh.ToPuzzleString();
        }

        public static string? ExportSquareAscii(Mesh mesh)
        {
            if (!TryGetSquareDimensions(mesh, out int width, out int height))
                return null;

            var sb = new StringBuilder();
            for (int j = 0; j <= height; j++)
            {
                // Horizontal line: intersections (O/*/+) and horizontal edges (-/x/ )
                for (int i = 0; i <= width; i++)
                {
                    int interIdx = i * (height + 1) + j;
                    sb.Append(mesh.Intersections[interIdx].Type switch
                    {
                        IntersType.Twist => 'O',
                        IntersType.Turn => '*',
                        _ => '+'
                    });
                    if (i < width)
                    {
                        int end = (i + 1) * (height + 1) + j;
                        EdgeState st = TryFindEdgeJoining(mesh, interIdx, end, out int edgeIdx)
                            ? mesh.Edges[edgeIdx].State
                            : EdgeState.Empty;
                        sb.Append(st == EdgeState.Filled ? '-' : (st == EdgeState.Excluded ? 'x' : ' '));
                    }
                }
                sb.AppendLine();

                // Vertical line: vertical edges (|/x/ ) separated by cell interiors (' ')
                if (j < height)
                {
                    for (int i = 0; i <= width; i++)
                    {
                        int start = i * (height + 1) + j;
                        int end = start + 1;
                        EdgeState st = TryFindEdgeJoining(mesh, start, end, out int edgeIdx)
                            ? mesh.Edges[edgeIdx].State
                            : EdgeState.Empty;
                        sb.Append(st == EdgeState.Filled ? '|' : (st == EdgeState.Excluded ? 'x' : ' '));
                        if (i < width)
                            sb.Append(' ');
                    }
                    sb.AppendLine();
                }
            }
            return sb.ToString();
        }

        public static bool TryParsePuzzle(string rawInput, out Mesh? mesh, out string detectedFormat, out int width, out int height)
        {
            mesh = null;
            detectedFormat = string.Empty;
            width = 0;
            height = 0;
            if (string.IsNullOrWhiteSpace(rawInput))
                return false;

            string trimmed = rawInput.Trim();

            // 1. Compact URL / #p= / #puzzle= hash / Pearl descriptor (wxh:...) / compact codec
            string codeCandidate = ExtractShareCode(trimmed);
            if (!trimmed.Contains('\n') && PuzzleCodec.TryDecode(codeCandidate, out var decodedMesh, out int dw, out int dh, out _) && decodedMesh != null)
            {
                mesh = decodedMesh;
                detectedFormat = (!codeCandidate.Contains('/') && codeCandidate.Contains(':'))
                    ? "Pearl Descriptor"
                    : "Share Link";
                width = dw;
                height = dh;
                return true;
            }

            string[] lines = rawInput.Split(LineSeparators, StringSplitOptions.RemoveEmptyEntries);
            if (lines.Length == 0)
                return false;

            // 2. Square Masyu ASCII grid (+, O, *, -, |, x)
            if (IsAsciiHeader(lines[0]) && TryParseSquareAscii(lines, out mesh, out width, out height))
            {
                detectedFormat = "ASCII Grid";
                return true;
            }

            // 3. Native .twist serialization
            try
            {
                var nativeMesh = new Mesh(0, 0, MeshType.Square);
                if (nativeMesh.LoadFromText(lines))
                {
                    mesh = nativeMesh;
                    detectedFormat = "Twist File";
                    if (TryGetSquareDimensions(nativeMesh, out int sw, out int sh))
                    {
                        width = sw;
                        height = sh;
                    }
                    return true;
                }
            }
            catch
            {
            }

            // 4. Legacy compressed share token fallback
            if (!trimmed.Contains('\n'))
            {
                try
                {
                    string decompressed = PuzzleCompression.Decompress(codeCandidate);
                    if (!string.IsNullOrWhiteSpace(decompressed))
                    {
                        var legacyMesh = new Mesh(0, 0, MeshType.Square);
                        if (legacyMesh.LoadFromText(decompressed.AsSpan()))
                        {
                            mesh = legacyMesh;
                            detectedFormat = "Legacy Share Link";
                            if (TryGetSquareDimensions(legacyMesh, out int lw, out int lh))
                            {
                                width = lw;
                                height = lh;
                            }
                            return true;
                        }
                    }
                }
                catch
                {
                }
            }

            return false;
        }

        private static bool IsAsciiHeader(string firstLine)
            => firstLine.Length >= 3 && firstLine[0] is '+' or 'O' or 'o' or '*' or '@';

        private static string ExtractShareCode(string trimmed)
        {
            int hashIdx = trimmed.IndexOf('#');
            string fragment = hashIdx >= 0 ? trimmed[(hashIdx + 1)..] : trimmed;
            if (fragment.StartsWith("p=", StringComparison.Ordinal))
                return fragment[2..];
            if (fragment.StartsWith("puzzle=", StringComparison.Ordinal))
                return fragment[7..];
            return fragment;
        }

        private static bool TryParseSquareAscii(string[] lines, out Mesh? mesh, out int width, out int height)
        {
            mesh = null;
            width = 0;
            height = 0;

            if (lines.Length < 3 || lines.Length % 2 == 0)
                return false;

            int lineLen = lines[0].Length;
            if (lineLen < 3 || lineLen % 2 == 0)
                return false;

            width = lineLen / 2;
            height = lines.Length / 2;
            if (width <= 0 || height <= 0 || width > 100 || height > 100)
                return false;

            var candidate = new Mesh(width, height, MeshType.Square);

            for (int i = 0; i < lines.Length; i++)
            {
                string line = lines[i];
                if (line.Length < lineLen)
                    line = line.PadRight(lineLen);

                int row = i / 2;
                if (i % 2 == 0)
                {
                    // Intersections at even positions
                    for (int j = 0; j < lineLen; j += 2)
                    {
                        int col = j / 2;
                        int interIdx = col * (height + 1) + row;
                        candidate.Intersections[interIdx].Type = line[j] switch
                        {
                            'O' or 'o' or 'W' or 'w' => IntersType.Twist,
                            '*' or '@' or 'B' or 'b' => IntersType.Turn,
                            _ => IntersType.Unknown
                        };
                    }

                    // Horizontal edges at odd positions
                    for (int j = 1; j < lineLen; j += 2)
                    {
                        char ch = line[j];
                        if (ch == ' ') continue;
                        int start = (j / 2) * (height + 1) + row;
                        int end = (j / 2 + 1) * (height + 1) + row;
                        if (TryFindEdgeJoining(candidate, start, end, out int edgeIdx))
                        {
                            EdgeState st = ch == '-' ? EdgeState.Filled : EdgeState.Excluded;
                            new SetAction(candidate, edgeIdx, st).Perform();
                        }
                    }
                }
                else
                {
                    // Vertical edges at even positions
                    for (int j = 0; j < lineLen; j += 2)
                    {
                        char ch = line[j];
                        if (ch == ' ') continue;
                        int start = (j / 2) * (height + 1) + row;
                        int end = start + 1;
                        if (TryFindEdgeJoining(candidate, start, end, out int edgeIdx))
                        {
                            EdgeState st = ch == '|' ? EdgeState.Filled : EdgeState.Excluded;
                            new SetAction(candidate, edgeIdx, st).Perform();
                        }
                    }
                }
            }

            mesh = candidate;
            return true;
        }
    }
}
