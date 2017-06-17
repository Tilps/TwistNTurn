using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bridge.Html5;
using TwistNTurn;

namespace TwistNTurnBridge
{
    class EdgeDisplay
    {

        public EdgeDisplay(CanvasRenderingContext2D context, double x1, double x2, double y1, double y2)
        {
            this.context = context;
            this.x1 = x1;
            this.x2 = x2;
            this.y1 = y1;
            this.y2 = y2;
        }

        private CanvasRenderingContext2D context;
        private double x1, x2, y1, y2;
        private bool dirty = true;



        public bool Marked
        {
            get { return marked; }
            set
            {
                if (value != marked)
                {
                    marked = value;
                    dirty = true;
                }
            }
        }

        private bool marked;

        public EdgeState EdgeState
        {
            get { return edgeState; }
            set
            {
                if (value != edgeState)
                {
                    edgeState = value;
                    dirty = true;
                }
            }
        }

        private EdgeState edgeState;

        public void Render()
        {
            if (!dirty) return;
            context.GlobalCompositeOperation = CanvasTypes.CanvasCompositeOperationType.DestinationOut;
            context.StrokeStyle = "black";
            context.LineWidth = 4;
            context.BeginPath();
            context.MoveTo(x1, y1);
            context.LineTo(x2, y2);
            context.Stroke();
            double midx = (x1 + x2) / 2.0;
            double midy = (y1 + y2) / 2.0;
            double len = Math.Sqrt(Math.Pow(x2 - x1, 2) + Math.Pow(y2 - y1, 2)) / 8;
            context.BeginPath();
            context.MoveTo(midx + len + 2, midy + len + 2);
            context.LineTo(midx - len - 2, midy - len - 2);
            context.MoveTo(midx + len + 2, midy - len - 2);
            context.LineTo(midx - len - 2, midy + len + 2);
            context.Stroke();

            context.GlobalCompositeOperation = CanvasTypes.CanvasCompositeOperationType.SourceOver;
            context.StrokeStyle = EdgeState == EdgeState.Filled ? (Marked ? "green" : "black") : "lightgrey";
            context.LineWidth = EdgeState == EdgeState.Filled ? 2 : 1;
            context.BeginPath();
            context.MoveTo(x1, y1);
            context.LineTo(x2, y2);
            context.Stroke();
            if (EdgeState == EdgeState.Excluded)
            {
                context.StrokeStyle = Marked ? "green" : "black";
                context.BeginPath();
                context.MoveTo(midx + len, midy + len);
                context.LineTo(midx - len, midy - len);
                context.MoveTo(midx + len, midy - len);
                context.LineTo(midx - len, midy + len);
                context.Stroke();
            }

        }

    }

}
