
    using System;

namespace Bridge.Html5
{

   /// <summary>
    /// External binding to the Worker HTML5 interface
    /// </summary>
    [External]
    [Name("Worker")]
    public
    class Worker
	{
		public Worker(string Uri)
		{
		}

        [Name("postMessage")]
		public extern void PostMessage(object o);

	    [Name("terminate")]
		public extern void Terminate();

        public class DataEvent
		{
		    [Name("data")]
			public object Data;
		}

		[Name("onmessage")]
		public Action<DataEvent> OnMessage;
	}
}
