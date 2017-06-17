using Bridge;
using Bridge.Html5;
using System;
using TwistNTurn;
using Newtonsoft.Json;

namespace TwistNTurnBridge
{
    public class App
    {
        private static LoopDisplay display;
        private static HTMLInputElement sizeInput;
        private static HTMLSelectElement difficultyPicker;
        private static Text timer;
        private static DateTime start = DateTime.MinValue;
        private static HTMLDivElement displayHost;
        private static HTMLButtonElement generateButton;
        public static void Main()
        {
            if (WorkerSpawn.CheckIsWorker())
            {
                return;
            }
            var meta = Document.CreateElement<HTMLMetaElement>("meta");
            meta.Name = "viewport";
            meta.Content = "width=device-width, initial-scale=1";
            Document.Head.AppendChild(meta);
            displayHost = Document.CreateElement<HTMLDivElement>("div");
            displayHost.Style.Height = "100%";
            displayHost.Style.Width = "100vw";
            displayHost.Style.Display = Display.TableRow;
            display = new LoopDisplay(displayHost);
            Document.Body.AppendChild(Document.CreateTextNode("Size:"));
            sizeInput = Document.CreateElement<HTMLInputElement>("input");
            sizeInput.Placeholder = "10x10";
            sizeInput.Type = InputType.Text;
            Document.Body.AppendChild(sizeInput);
            difficultyPicker = Document.CreateElement<HTMLSelectElement>("select");
            difficultyPicker.Add(Document.CreateElement<HTMLOptionElement>("option"));
            difficultyPicker.Add(Document.CreateElement<HTMLOptionElement>("option"));
            difficultyPicker.Add(Document.CreateElement<HTMLOptionElement>("option"));
            difficultyPicker.Add(Document.CreateElement<HTMLOptionElement>("option"));
            difficultyPicker.Add(Document.CreateElement<HTMLOptionElement>("option"));
            difficultyPicker.Add(Document.CreateElement<HTMLOptionElement>("option"));
            difficultyPicker.Options[0].Text = "Trivial";
            difficultyPicker.Options[1].Text = "Easy";
            difficultyPicker.Options[2].Text = "Easyish";
            difficultyPicker.Options[3].Text = "Moderate";
            difficultyPicker.Options[4].Text = "Maybe Harder";
            difficultyPicker.Options[5].Text = "Unlimited";
            difficultyPicker.SelectedIndex = 1;
            Document.Body.AppendChild(difficultyPicker);
            generateButton = Document.CreateElement<HTMLButtonElement>("button");
            generateButton.TextContent = "Generate";
            generateButton.OnClick = OnClick;
            Document.Body.AppendChild(generateButton);
            timer = Document.CreateTextNode("0:00");
            Document.Body.AppendChild(timer);
            Document.Body.AppendChild(Document.CreateElement<HTMLBRElement>("br"));
            var undoButton = Document.CreateElement<HTMLButtonElement>("button");
            undoButton.TextContent = "Undo";
            undoButton.OnClick = e => display.Undo();
            Document.Body.AppendChild(undoButton);
            var redoButton = Document.CreateElement<HTMLButtonElement>("button");
            redoButton.TextContent = "Redo";
            redoButton.OnClick = e => display.Redo();
            Document.Body.AppendChild(redoButton);
            Document.Body.AppendChild(Document.CreateElement<HTMLBRElement>("br"));
            var fixButton = Document.CreateElement<HTMLButtonElement>("button");
            fixButton.TextContent = "Fix";
            fixButton.OnClick = e => display.Fix();
            Document.Body.AppendChild(fixButton);
            var unfixButton = Document.CreateElement<HTMLButtonElement>("button");
            unfixButton.TextContent = "Unfix";
            unfixButton.OnClick = e => display.Unfix();
            Document.Body.AppendChild(unfixButton);
            var revertButton = Document.CreateElement<HTMLButtonElement>("button");
            revertButton.TextContent = "Revert";
            revertButton.OnClick = e => display.RevertToFix();
            Document.Body.AppendChild(revertButton);



            Document.Body.AppendChild(displayHost);
            Document.Body.Style.Height = "100vh";
            Document.Body.Style.Width = "100vw";
            Document.Body.Style.Margin = "0";
            Document.Body.Style.Padding = "0";
            Document.Body.Style.Display = Display.Table;
            Document.DocumentElement.Style.Height = "100vh";
            Document.DocumentElement.Style.Width = "100vw";
            Document.DocumentElement.Style.Margin = "0";
            Document.DocumentElement.Style.Padding = "0";
            display.Mesh = new Mesh(10, 10, MeshType.Square);
            //display.Mesh.Generate();
            display.Mesh = display.Mesh;
            Window.OnResize = App_Resize;
            Window.OnKeyDown = display.OnKeyDown;
            Document.Body.OnLoad = App_Resize;
            Window.SetTimeout(App_Tick, 100);
        }

        public static void App_Tick()
        {
            if (start != DateTime.MinValue)
            {
                TimeSpan duration = DateTime.UtcNow - start;
                timer.TextContent = duration.ToString();
            }
            Window.SetTimeout(App_Tick, 100);
        }

        public static MeshType MeshTypeFromString(string typeName)
        {
            if (string.IsNullOrEmpty(typeName))
                return MeshType.Square;
            MeshType type = MeshType.Square;
            if (typeName == "Square")
                type = MeshType.Square;
            else if (typeName == "Square Symmetrical")
                type = MeshType.SquareSymmetrical;
            else if (typeName == "Triangle")
                type = MeshType.Triangle;
            else if (typeName == "Hexagon")
                type = MeshType.Hexagonal;
            else if (typeName == "Hexagon2")
                type = MeshType.Hexagonal2;
            else if (typeName == "Hexagon3")
                type = MeshType.Hexagonal3;
            else if (typeName == "Octagon")
                type = MeshType.Octagon;
            else if (typeName == "Square2")
                type = MeshType.Square2;
            else if (typeName == "Pentagon")
                type = MeshType.Pentagon;
            return type;
        }

        private static Worker worker;
        private static HTMLProgressElement progressBar;
        private static void OnClick(MouseEvent<HTMLButtonElement> mouseEvent)
        {
            MeshType type = MeshType.Square;
            int width;
            int height;
            if (!ParseSize(sizeInput.Value, type, out width, out height))
                return;
            Mesh mesh = MakeMesh(width, height, type, difficultyPicker.SelectedIndex);
            display.Mesh = mesh;
            
            worker = new Worker(Extensions.Window().URL.createObjectURL(new Blob(new BlobDataObject[]
            {
                @"
self.onmessage = function(e) { 
  if (e.data.href) {
    try { 
      importScripts(e.data.href);
    } catch (error) {
      console.log(e.data.href);  
      console.log(error);
    }
  } else {
    TwistNTurnBridge.WorkerSpawn.WorkerSpawn_OnMessage(e);
  }
}"
            }, new BlobPropertyBag() {Type="text/javascript"})));
            progressBar = Document.CreateElement<HTMLProgressElement>("progress");
            progressBar.Max = mesh.Intersections.Count;
            progressBar.Style.Position = Position.Absolute;
            progressBar.Style.Margin = "auto";
            progressBar.Style.Top = "0";
            progressBar.Style.Bottom = "0";
            progressBar.Style.Left = "0";
            progressBar.Style.Right = "0";
            progressBar.Style.ZIndex = "100";
            progressBar.Style.BackgroundColor = "white";
            displayHost.AppendChild(progressBar);
            generateButton.Disabled = true;
            worker.OnMessage += AppWorker_OnMessage;
            string to_load = Window.Location.Href.Substring(0, Window.Location.Href.LastIndexOf('/') + 1) + "bridge.min.js";
            worker.PostMessage(new InitialReuqest() { href = to_load});
            to_load = Window.Location.Href.Substring(0, Window.Location.Href.LastIndexOf('/') + 1) + "bridge.console.min.js";
            worker.PostMessage(new InitialReuqest() { href = to_load });
            to_load = Window.Location.Href.Substring(0, Window.Location.Href.LastIndexOf('/') + 1) + "bridge.meta.min.js";
            worker.PostMessage(new InitialReuqest() { href = to_load });
            to_load = Window.Location.Href.Substring(0, Window.Location.Href.LastIndexOf('/') + 1) + "newtonsoft.json.min.js";
            worker.PostMessage(new InitialReuqest() { href = to_load });
            to_load = Window.Location.Href.Substring(0, Window.Location.Href.LastIndexOf('/') + 1) + "TwistNTurnBridge.min.js";
            worker.PostMessage(new InitialReuqest() { href = to_load });
            to_load = Window.Location.Href.Substring(0, Window.Location.Href.LastIndexOf('/') + 1) + "TwistNTurnBridge.meta.min.js";
            worker.PostMessage(new InitialReuqest() { href = to_load });
            worker.PostMessage(JsonConvert.SerializeObject(new GenerateRequest() {Width=width, Height=height, Type=type, Difficulty = difficultyPicker.SelectedIndex}));
        }

        private static void AppWorker_OnMessage(Worker.DataEvent dataEvent)
        {
            string msg = (string) dataEvent.Data;
            if (msg.Length < 5)
            {
                progressBar.Value = int.Parse(msg);
            }
            else
            {
                Mesh current = display.Mesh;
                Mesh newMesh = new Mesh(0, 0, current.MeshType);
                newMesh.LoadFromText(msg.Split(new char[] {'\n', '\r'}, StringSplitOptions.RemoveEmptyEntries));
                display.Mesh = newMesh;
                progressBar.Remove();
                progressBar = null;
                generateButton.Disabled = false;
                start = DateTime.UtcNow;
            }
        }

        public static Mesh MakeMesh(int width, int height, MeshType type, int difficulty)
        {
            Mesh mesh = new Mesh(width, height, type);
            mesh.ConsiderMultipleLoops = difficulty > 0;
            mesh.IterativeRecMaxDepth = 1;
            if (difficulty > 1)
            {
                mesh.UseCellColoring = true;
                mesh.UseCellColoringTrials = true;
                mesh.UseColoring = true;
                mesh.UseEdgeRestricts = true;
                mesh.UseDerivedColoring = true;
                mesh.UseMerging = true;
                mesh.UseCellPairsTopLevel = true;
            }
            if (difficulty < 5)
            {
                mesh.IterativeSolverDepth = Math.Max(0, difficulty - 2);
            }
            else
            {
                mesh.SolverMethod = SolverMethod.Recursive;
            }
            mesh.GenerateBoringFraction = 0.01;
            return mesh;
        }

        public static bool ParseSize(string val, MeshType type, out int width, out int height)
        {
            val = val.Trim();
            width = 10;
            height = 10;
            if (val.Length == 0)
            {
                if (type == MeshType.Octagon)
                {
                    width = 5;
                    height = 5;
                }
                else if (type == MeshType.Square2)
                {
                    width = 5;
                    height = 5;
                }
                else if (type == MeshType.Hexagonal)
                {
                    width = 5;
                    height = 10;
                }
                else if (type == MeshType.Hexagonal2)
                {
                    width = 6;
                    height = 6;
                }
                else if (type == MeshType.Hexagonal3)
                {
                    width = 4;
                    height = 4;
                }
                else if (type == MeshType.Triangle)
                {
                    width = 6;
                    height = 6;
                }
                else if (type == MeshType.Pentagon)
                {
                    width = 6;
                    height = 6;
                }
                return true;
            }
            string[] bits = val.Split('x');
            if (bits.Length == 2)
            {
                if (!int.TryParse(bits[0], out width))
                    return false;
                if (!int.TryParse(bits[1], out height))
                    return false;
                return true;
            }
            else if (bits.Length == 1)
            {
                if (!int.TryParse(bits[0], out width))
                    return false;
                height = width;
                return true;
            }
            return false;

        }

        private static void App_Resize(Event arg)
        {
            display.OnScaleChanged();
        }

    }

    [Reflectable]
    public class GenerateRequest
    {
        public int Width;
        public int Height;
        public MeshType Type;
        public int Difficulty;
    }

    public class InitialReuqest
    {
        public string href;
    }


}