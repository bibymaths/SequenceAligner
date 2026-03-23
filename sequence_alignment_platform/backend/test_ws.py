import asyncio
import httpx
import websockets

async def run_full_test():
    base_url = "http://127.0.0.1:8000"
    ws_url = "ws://127.0.0.1:8000/ws/logs"

    # 1. Trigger the Alignment (POST)
    print("🚀 Sending files to /align...")
    async with httpx.AsyncClient() as client:
        # Create dummy files in memory for the test
        files = {
            "query": ("query.fa", ">seq1\nACGT", "text/plain"),
            "target": ("target.fa", ">seq2\nACGT", "text/plain")
        }
        response = await client.post(f"{base_url}/align", files=files)

        if response.status_code != 200:
            print(f"❌ Failed to start: {response.text}")
            return

        session_id = response.json()["session_id"]
        print(f"✅ Session Created: {session_id}")

    # 2. Immediately Connect to WebSocket
    print(f"🔗 Connecting to WebSocket for logs...")
    async with websockets.connect(f"{ws_url}/{session_id}") as ws:
        try:
            # We use a timeout so the script doesn't hang forever
            while True:
                # Wait for a message, but if nothing comes for 5s, assume we're done
                message = await asyncio.wait_for(ws.recv(), timeout=5.0)
                print(f"📥 {message.strip()}")
                if "Analysis complete" in message or "completed" in message:
                    print("🏁 Alignment finished successfully!")
                    break
        except asyncio.TimeoutError:
            print("🕒 No more logs received (Timeout).")

if __name__ == "__main__":
    asyncio.run(run_full_test())