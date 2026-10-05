using System;
using System.IO;
using System.IO.Compression;
using System.Text;

namespace TwistNTurn
{
    public static class PuzzleCompression
    {
        public static string Compress(string text) => CompressToUrlSafe(text);
        public static string Decompress(string urlSafeBase64) => DecompressFromUrlSafe(urlSafeBase64);

        public static string CompressToUrlSafe(string text)
        {
            if (string.IsNullOrEmpty(text)) return string.Empty;
            byte[] bytes = Encoding.UTF8.GetBytes(text);
            using var ms = new MemoryStream();
            using (var def = new DeflateStream(ms, CompressionLevel.Optimal, true))
            {
                def.Write(bytes, 0, bytes.Length);
            }
            return Convert.ToBase64String(ms.ToArray())
                .Replace('+', '-')
                .Replace('/', '_')
                .TrimEnd('=');
        }

        public static string DecompressFromUrlSafe(string urlSafeBase64)
        {
            if (string.IsNullOrEmpty(urlSafeBase64)) return string.Empty;
            string padded = urlSafeBase64.Replace('-', '+').Replace('_', '/');
            switch (padded.Length % 4)
            {
                case 2: padded += "=="; break;
                case 3: padded += "="; break;
            }
            byte[] compressed = Convert.FromBase64String(padded);
            using var inMs = new MemoryStream(compressed);
            using var def = new DeflateStream(inMs, CompressionMode.Decompress);
            using var outMs = new MemoryStream();
            def.CopyTo(outMs);
            return Encoding.UTF8.GetString(outMs.ToArray());
        }
    }
}
