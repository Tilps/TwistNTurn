using System.Text.Json.Serialization;

namespace TwistNTurn.Web.Services
{
    [JsonSerializable(typeof(SavedGameState))]
    [JsonSerializable(typeof(UserPreferences))]
    internal partial class TwistJsonContext : JsonSerializerContext
    {
    }
}
