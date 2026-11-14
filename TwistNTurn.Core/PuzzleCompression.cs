using System;
using System.Buffers.Text;
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
            using (var def = new DeflateStream(ms, CompressionLevel.Optimal, leaveOpen: true))
            {
                def.Write(bytes);
            }
            return ms.TryGetBuffer(out var buffer)
                ? Base64Url.EncodeToString(buffer)
                : Base64Url.EncodeToString(ms.ToArray());
        }

        public static string DecompressFromUrlSafe(string urlSafeBase64)
        {
            if (string.IsNullOrEmpty(urlSafeBase64)) return string.Empty;
            byte[] compressed = Base64Url.DecodeFromChars(urlSafeBase64.AsSpan().TrimEnd('='));
            using var inMs = new MemoryStream(compressed, writable: false);
            using var def = new DeflateStream(inMs, CompressionMode.Decompress);
            using var outMs = new MemoryStream();
            def.CopyTo(outMs);
            return outMs.TryGetBuffer(out var buffer)
                ? Encoding.UTF8.GetString(buffer)
                : Encoding.UTF8.GetString(outMs.ToArray());
        }
    }
}
