using System;
using System.Buffers.Text;
using System.IO;
using System.IO.Compression;
using System.Text;

namespace TwistNTurn
{
    public static class PuzzleCompression
    {
        public const int MaxDecompressedBytes = 512 * 1024;

        public static string Compress(string text) => CompressToUrlSafe(text);
        public static string Decompress(string urlSafeBase64) => DecompressFromUrlSafe(urlSafeBase64);

        public static string CompressToUrlSafe(string text)
        {
            if (string.IsNullOrEmpty(text)) return string.Empty;
            using var input = new StringStream(text, Encoding.UTF8);
            using var ms = new MemoryStream();
            using (var def = new DeflateStream(ms, CompressionLevel.Optimal, leaveOpen: true))
            {
                input.CopyTo(def);
            }
            ReadOnlySpan<byte> compressed = ms.TryGetBuffer(out var buffer) ? buffer : ms.ToArray();
            return Base64Url.EncodeToString(compressed);
        }

        public static bool TryDecompressFromUrlSafe(string urlSafeBase64, out string result, int maxBytes = MaxDecompressedBytes)
        {
            result = string.Empty;
            if (string.IsNullOrEmpty(urlSafeBase64))
                return true;
            if (maxBytes <= 0 || urlSafeBase64.Length > maxBytes * 2)
                return false;

            try
            {
                byte[] compressed = Base64Url.DecodeFromChars(urlSafeBase64.AsSpan().TrimEnd('='));
                using var inMs = new ReadOnlyMemoryStream(compressed);
                using var def = new DeflateStream(inMs, CompressionMode.Decompress);
                using var outMs = new MemoryStream();
                byte[] chunk = new byte[4096];
                int totalRead = 0;
                int read;
                while ((read = def.Read(chunk, 0, chunk.Length)) > 0)
                {
                    totalRead += read;
                    if (totalRead > maxBytes)
                        return false;
                    outMs.Write(chunk, 0, read);
                }
                ReadOnlySpan<byte> utf8 = outMs.TryGetBuffer(out var buffer) ? buffer : outMs.ToArray();
                result = Encoding.UTF8.GetString(utf8);
                return true;
            }
            catch (Exception ex) when (ex is FormatException or InvalidDataException or ArgumentException)
            {
                return false;
            }
        }

        public static string DecompressFromUrlSafe(string urlSafeBase64)
        {
            if (string.IsNullOrEmpty(urlSafeBase64)) return string.Empty;
            if (!TryDecompressFromUrlSafe(urlSafeBase64, out string result))
                throw new InvalidDataException("Compressed puzzle payload is invalid or exceeds maximum allowed size.");
            return result;
        }
    }
}
