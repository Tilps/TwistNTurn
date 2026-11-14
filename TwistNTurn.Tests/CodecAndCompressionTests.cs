using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using TwistNTurn;

namespace TwistNTurn.Tests
{
    [TestClass]
    public class CodecAndCompressionTests
    {
        [TestMethod]
        [DataRow(MeshType.Square, 5, 5)]
        [DataRow(MeshType.SquareSymmetrical, 6, 6)]
        [DataRow(MeshType.Hexagonal2, 5, 5)]
        [DataRow(MeshType.Hexagonal3, 3, 3)]
        [DataRow(MeshType.Hexagonal4, 3, 3)]
        public void PuzzleCodec_RoundTrip_PreservesCluesAndTopology(MeshType type, int width, int height)
        {
            var original = new Mesh(width, height, type);

            // Set a few arbitrary clues
            if (original.Intersections.Count > 4)
            {
                original.Intersections[1].Type = IntersType.Twist; // white pearl
                original.Intersections[3].Type = IntersType.Turn;  // black pearl
            }

            string encoded = PuzzleCodec.Encode(original, width, height);
            Assert.IsFalse(string.IsNullOrWhiteSpace(encoded));

            var decoded = PuzzleCodec.Decode(encoded);
            Assert.IsNotNull(decoded);
            Assert.AreEqual(original.MeshType, decoded.MeshType);
            Assert.AreEqual(original.Intersections.Count, decoded.Intersections.Count);

            for (int i = 0; i < original.Intersections.Count; i++)
            {
                Assert.AreEqual(original.Intersections[i].Type, decoded.Intersections[i].Type,
                    $"Clue mismatch at intersection {i}");
            }
        }

        [TestMethod]
        public void PuzzleCodec_LegacyFormat_DecodesCorrectly()
        {
            // Legacy format: "{width}x{height}:{cluestring}"
            // Clue string contains '.' for unknown, '1' for twist (white), '2' for turn (black)
            string legacy = "3x3:.1.2.....";
            var mesh = PuzzleCodec.Decode(legacy);

            Assert.IsNotNull(mesh);
            Assert.AreEqual(MeshType.Square, mesh.MeshType);
            Assert.AreEqual(IntersType.Twist, mesh.Intersections[1].Type);
            Assert.AreEqual(IntersType.Turn, mesh.Intersections[3].Type);
            Assert.AreEqual(IntersType.Unknown, mesh.Intersections[0].Type);
        }

        [TestMethod]
        [DataRow("a")]
        [DataRow("ab")]
        [DataRow("abc")]
        [DataRow("abcd")]
        [DataRow("sq/10x10/5w3b12w8b...sample_puzzle_string_with_data")]
        public void PuzzleCompression_RoundTrip_CompressesAndDecompresses(string original)
        {
            string compressed = PuzzleCompression.Compress(original);

            Assert.IsFalse(string.IsNullOrWhiteSpace(compressed));
            Assert.IsFalse(compressed.Contains('+'));
            Assert.IsFalse(compressed.Contains('/'));
            Assert.IsFalse(compressed.Contains('='));
            string decompressed = PuzzleCompression.Decompress(compressed);

            Assert.AreEqual(original, decompressed);
        }

        [TestMethod]
        public void EdgeStateSerialization_RoundTrip_PreservesAllStates()
        {
            var mesh = new Mesh(4, 4, MeshType.Square);

            mesh.Edges[0].State = EdgeState.Filled;
            mesh.Edges[1].State = EdgeState.Excluded;
            mesh.Edges[2].State = EdgeState.Empty;
            mesh.Edges[3].State = EdgeState.Filled;

            string stateString = mesh.ToEdgeStateString();
            Assert.IsFalse(string.IsNullOrWhiteSpace(stateString));

            var mesh2 = new Mesh(4, 4, MeshType.Square);
            mesh2.FromEdgeStateString(stateString);

            for (int i = 0; i < mesh.Edges.Count; i++)
            {
                Assert.AreEqual(mesh.Edges[i].State, mesh2.Edges[i].State,
                    $"Edge state mismatch at edge {i}");
            }
        }

        [TestMethod]
        [DataRow("Hexagon 2", MeshType.Hexagonal2)]
        [DataRow("Hexagon2", MeshType.Hexagonal2)]
        [DataRow("Hexagonal2", MeshType.Hexagonal2)]
        [DataRow("hex2", MeshType.Hexagonal2)]
        [DataRow("Hexagon 3", MeshType.Hexagonal3)]
        [DataRow("Hexagon3", MeshType.Hexagonal3)]
        [DataRow("Hexagonal3", MeshType.Hexagonal3)]
        [DataRow("hex3", MeshType.Hexagonal3)]
        [DataRow("Hexagon 4", MeshType.Hexagonal4)]
        [DataRow("Hexagon4", MeshType.Hexagonal4)]
        [DataRow("Hexagonal4", MeshType.Hexagonal4)]
        [DataRow("hex4", MeshType.Hexagonal4)]
        [DataRow("Square", MeshType.Square)]
        [DataRow("sq", MeshType.Square)]
        [DataRow("Square Symmetrical", MeshType.SquareSymmetrical)]
        [DataRow("SquareSymmetrical", MeshType.SquareSymmetrical)]
        [DataRow("sqs", MeshType.SquareSymmetrical)]
        public void PuzzleHelper_MeshTypeFromString_ParsesCorrectly(string input, MeshType expected)
        {
            var result = PuzzleHelper.MeshTypeFromString(input);
            Assert.AreEqual(expected, result);
        }

        [TestMethod]
        [DataRow(MeshType.Hexagonal2, 5, 5)]
        [DataRow(MeshType.Hexagonal3, 3, 3)]
        [DataRow(MeshType.Hexagonal4, 3, 3)]
        public void PuzzleHelper_HexagonalTopologies_DefaultSizeAndParse(MeshType type, int expectedW, int expectedH)
        {
            PuzzleHelper.GetDefaultSize(type, out int w, out int h);
            Assert.AreEqual(expectedW, w);
            Assert.AreEqual(expectedH, h);

            bool ok = PuzzleHelper.ParseSize($"{expectedW}x{expectedH}", type, out int pw, out int ph);
            Assert.IsTrue(ok);
            Assert.AreEqual(expectedW, pw);
            Assert.AreEqual(expectedH, ph);
        }
    }
}

