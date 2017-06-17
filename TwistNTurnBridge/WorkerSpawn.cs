using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bridge;
using Bridge.Html5;
using TwistNTurn;
using Newtonsoft.Json;

namespace TwistNTurnBridge
{
    class WorkerSpawn
    {
        public static void WorkerSpawn_OnMessage(Worker.DataEvent dataEvent)
        {
            GenerateRequest request = JsonConvert.DeserializeObject<GenerateRequest>((string) dataEvent.Data);
            Mesh m = App.MakeMesh(request.Width, request.Height, request.Type, request.Difficulty);
            m.PrunedCountProgress += WorkerSpawn_OnPrunedCountProgress;
            counter = 0;
            m.Generate();
            
            StringBuilder builder = new StringBuilder();
            m.Save(builder);
            PostMessage(builder.ToString());
        }

        private static int counter = 0;

        private static void WorkerSpawn_OnPrunedCountProgress(object sender, EventArgs eventArgs)
        {
            counter++;
            PostMessage(counter.ToString());
        }

        [Template("postMessage({data})")]
        public static void PostMessage(object data)
        {
        }

        [Template("typeof window == 'undefined'")]
        public static bool CheckIsWorker()
        {
            return false;
        }

    }

    public class URLType
    {

        public extern string createObjectURL(Blob b);
    }

    public class WindowGlobal
    {
        public URLType URL {
            get { return null; }
        }
    }

    public class SelfGlobal
    {
        [Name("onmessage")]
        public Action<Worker.DataEvent> OnMessage;

    }

    public class ConsoleGlobal
    {
        [Name("log")]
        public extern void Log(object message);
    }


    public class Extensions
    {
        [Template("window")]
        public static WindowGlobal Window()
        {
            return null;
        }
        [Template("self")]
        public static SelfGlobal Self()
        {
            return null;
        }

        [Template("console")]
        public static ConsoleGlobal Console()
        {
            return null;
        }
    }
}
