using System.Text.Json.Serialization;

namespace TwistNTurn.Web.Services
{
    [JsonSerializable(typeof(SavedGameState))]
    [JsonSerializable(typeof(UserPreferences))]
    public partial class TwistJsonContext : JsonSerializerContext
    {
    }
}
