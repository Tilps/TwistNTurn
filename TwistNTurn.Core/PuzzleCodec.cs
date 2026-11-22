using System;
using System.Globalization;
using System.Text;

namespace TwistNTurn
{
    public static class PuzzleCodec
    {
        public static string ToTypeCode(MeshType type) => type switch
        {
            MeshType.Square => "sq",
            MeshType.SquareSymmetrical => "sqs",
            MeshType.Hexagonal2 => "hex2",
            MeshType.Hexagonal3 => "hex3",
            MeshType.Hexagonal4 => "hex4",
            _ => "sq"
        };

        public static bool TryFromTypeCode(string code, out MeshType type)
        {
            if (string.IsNullOrWhiteSpace(code))
            {
                type = MeshType.Square;
                return false;
            }

            switch (code.Trim().ToLowerInvariant())
            {
                case "sq":
                case "square":
                    type = MeshType.Square; return true;
                case "sqs":
                case "squaresymmetrical":
                    type = MeshType.SquareSymmetrical; return true;
                case "hex2":
                case "hexagonal2":
                case "hexagon2":
                    type = MeshType.Hexagonal2; return true;
                case "hex3":
                case "hexagonal3":
                case "hexagon3":
                    type = MeshType.Hexagonal3; return true;
                case "hex4":
                case "hexagonal4":
                case "hexagon4":
                    type = MeshType.Hexagonal4; return true;
                default:
                    type = MeshType.Square;
                    return false;
            }
        }

        public static MeshType FromTypeCode(string code) => TryFromTypeCode(code, out var type) ? type : MeshType.Square;

        public static string Encode(Mesh mesh, int width, int height)
        {
            if (mesh == null) return string.Empty;

            var sb = new StringBuilder();
            sb.Append(ToTypeCode(mesh.MeshType));
            sb.Append('/');
            sb.Append(width.ToString(CultureInfo.InvariantCulture));
            sb.Append('x');
            sb.Append(height.ToString(CultureInfo.InvariantCulture));
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
            if (code.StartsWith("p=", StringComparison.Ordinal)) code = code[2..];
            else if (code.StartsWith("puzzle=", StringComparison.Ordinal)) code = code[7..];

            // Legacy format: {w}x{h}:{cluestring}
            int colonIndex = code.IndexOf(':');
            if (colonIndex > 0 && !code.Contains('/'))
            {
                try
                {
                    string[] dimParts = code[..colonIndex].Split('x');
                    if (dimParts.Length != 2
                        || !int.TryParse(dimParts[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out int w)
                        || !int.TryParse(dimParts[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int h)
                        || w < 2 || h < 2
                        || w - 1 > PuzzleHelper.MaxGridDimension || h - 1 > PuzzleHelper.MaxGridDimension)
                    {
                        return false;
                    }

                    mesh = Mesh.FromPuzzleString(code);
                    width = w - 1;
                    height = h - 1;
                    type = MeshType.Square;
                    return true;
                }
                catch
                {
                    mesh = null;
                    return false;
                }
            }

            string[] parts = code.Split('/', StringSplitOptions.None);
            if (parts.Length is < 2 or > 3) return false;

            if (!TryFromTypeCode(parts[0], out type))
                return false;

            string sizeStr = parts[1];
            if (string.IsNullOrWhiteSpace(sizeStr) || !PuzzleHelper.ParseSize(sizeStr, type, out width, out height))
                return false;

            try
            {
                mesh = new Mesh(width, height, type);
                if (parts.Length == 3 && !string.IsNullOrEmpty(parts[2]))
                {
                    string clueString = parts[2];
                    int intersIndex = 0;
                    for (int i = 0; i < clueString.Length; i++)
                    {
                        char c = clueString[i];
                        if (c == '0')
                        {
                            if (intersIndex >= mesh.Intersections.Count)
                            {
                                mesh = null;
                                return false;
                            }
                            mesh.Intersections[intersIndex].Type = IntersType.Twist;
                            intersIndex++;
                        }
                        else if (c == '1')
                        {
                            if (intersIndex >= mesh.Intersections.Count)
                            {
                                mesh = null;
                                return false;
                            }
                            mesh.Intersections[intersIndex].Type = IntersType.Turn;
                            intersIndex++;
                        }
                        else if (c is >= 'a' and <= 'z')
                        {
                            intersIndex += (c - 'a' + 1);
                            if (intersIndex > mesh.Intersections.Count)
                            {
                                mesh = null;
                                return false;
                            }
                        }
                        else
                        {
                            mesh = null;
                            return false;
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
