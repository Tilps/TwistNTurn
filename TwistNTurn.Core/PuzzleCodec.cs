#nullable enable
using System;
using System.Text;

namespace TwistNTurn
{
    public static class PuzzleCodec
    {
        public static string ToTypeCode(MeshType type) => type switch
        {
            MeshType.Square => "sq",
            MeshType.SquareSymmetrical => "sqs",
            MeshType.Hexagonal3 => "hex3",
            _ => "sq"
        };

        public static MeshType FromTypeCode(string code) => code.ToLowerInvariant() switch
        {
            "sq" or "square" => MeshType.Square,
            "sqs" or "squaresymmetrical" => MeshType.SquareSymmetrical,
            "hex3" or "hexagonal3" or "hexagon3" => MeshType.Hexagonal3,
            _ => MeshType.Square
        };

        public static string Encode(Mesh mesh, int width, int height)
        {
            if (mesh == null) return string.Empty;

            var sb = new StringBuilder();
            sb.Append(ToTypeCode(mesh.MeshType));
            sb.Append('/');
            sb.Append(width);
            sb.Append('x');
            sb.Append(height);
            sb.Append('/');

            int emptyRun = 0;
            for (int i = 0; i < mesh.Intersections.Count; i++)
            {
                var type = mesh.Intersections[i].Type;
                if (type == IntersType.Unknown)
                {
                    emptyRun++;
                }
                else
                {
                    while (emptyRun > 26)
                    {
                        sb.Append('z');
                        emptyRun -= 26;
                    }
                    if (emptyRun > 0)
                    {
                        sb.Append((char)('a' + emptyRun - 1));
                        emptyRun = 0;
                    }

                    if (type == IntersType.Twist)
                        sb.Append('0');
                    else if (type == IntersType.Turn)
                        sb.Append('1');
                }
            }
            // Trailing empty intersections do not need to be written.

            return sb.ToString();
        }

        public static Mesh? Decode(string code)
        {
            if (TryDecode(code, out var mesh, out _, out _, out _))
                return mesh;
            return null;
        }

        public static bool TryDecode(string code, out Mesh? mesh, out int width, out int height, out MeshType type)
        {
            mesh = null;
            width = 0;
            height = 0;
            type = MeshType.Square;

            if (string.IsNullOrWhiteSpace(code)) return false;

            code = code.Trim();
            if (code.StartsWith('#')) code = code[1..];
            if (code.StartsWith("p=")) code = code[2..];
            else if (code.StartsWith("puzzle=")) code = code[7..];

            // Check if it's the legacy format: {w}x{h}:{cluestring}
            int colonIndex = code.IndexOf(':');
            if (colonIndex > 0 && !code.Contains('/'))
            {
                try
                {
                    mesh = Mesh.FromPuzzleString(code);
                    string[] dimParts = code.Substring(0, colonIndex).Split('x');
                    if (dimParts.Length == 2 && int.TryParse(dimParts[0], out int w) && int.TryParse(dimParts[1], out int h))
                    {
                        width = w - 1;
                        height = h - 1;
                    }
                    type = MeshType.Square;
                    return true;
                }
                catch
                {
                    return false;
                }
            }

            string[] parts = code.Split('/', StringSplitOptions.None);
            if (parts.Length < 2) return false;

            type = FromTypeCode(parts[0]);

            string sizeStr = parts[1];
            if (!PuzzleHelper.ParseSize(sizeStr, type, out width, out height))
            {
                return false;
            }

            try
            {
                mesh = new Mesh(width, height, type);
                if (parts.Length >= 3 && !string.IsNullOrEmpty(parts[2]))
                {
                    string clueString = parts[2];
                    int intersIndex = 0;
                    for (int i = 0; i < clueString.Length && intersIndex < mesh.Intersections.Count; i++)
                    {
                        char c = clueString[i];
                        if (c == '0')
                        {
                            mesh.Intersections[intersIndex].Type = IntersType.Twist;
                            intersIndex++;
                        }
                        else if (c == '1')
                        {
                            mesh.Intersections[intersIndex].Type = IntersType.Turn;
                            intersIndex++;
                        }
                        else if (c >= 'a' && c <= 'z')
                        {
                            intersIndex += (c - 'a' + 1);
                        }
                    }
                }
                return true;
            }
            catch
            {
                mesh = null;
                return false;
            }
        }
    }
}
