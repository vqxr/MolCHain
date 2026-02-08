# MolChain Dashboard Troubleshooting Guide

## The Dashboard Shows "Loading..."

This means the frontend can't connect to the backend. Follow these steps:

### Step 1: Start the Flask Backend

Open a terminal and run:

```bash
cd /Users/abdurrahman/Documents/MolChain/backend
python3 app.py
```

You should see output like:
```
🚀 Deploying MolChain Backend...
✅ Backend initialized successfully

📊 Dashboard Statistics:
  • Total Molecules: 2
  • Total Contributions: 2
  • AI Agents: 5
  ...

🎉 All Systems Initialized!

📡 Essential Endpoints:
  Health Check:      http://localhost:5000/api/health
  All Molecules:     http://localhost:5000/api/molecules
  AI Analysis:       POST http://localhost:5000/api/analyze
  ...

✅ Server ready! React frontend should connect automatically
⏰ Press Ctrl+C to stop the server
```

### Step 2: Verify Backend is Running

In another terminal, test the health endpoint:

```bash
curl http://localhost:5000/api/health
```

You should get JSON response:
```json
{
  "status": "healthy",
  "service": "MolChain AI + Blockchain Backend",
  "version": "2.0.0-HACKATHON"
}
```

### Step 3: Start the Frontend

In yet another terminal:

```bash
cd /Users/abdurrahman/Documents/MolChain/frontend
npm run dev
```

You should see:
```
  VITE v5... ready in XXX ms

  ➜  Local:   http://localhost:5173/
```

### Step 4: Open the Dashboard

Go to: `http://localhost:5173`

---

## Still Seeing "Loading..."?

### Check 1: Is the Backend Running?
```bash
# In a terminal, check if port 5000 is in use
lsof -i :5000
# OR
netstat -an | grep 5000
```

❌ No output? → Backend is NOT running  
✅ Shows `python` or `python3` → Backend IS running

### Check 2: Is the Frontend Connecting?
Open your browser's **Developer Tools** (F12 or Cmd+Option+I):

1. Go to **Console** tab
2. Look for errors like:
   - `Failed to connect to http://localhost:5000...`
   - `CORS error`
   - `Connection refused`

### Check 3: CORS Issues?
If you see CORS errors, the backend needs to allow the frontend origin.

Edit `backend/app.py` line ~24:

```python
CORS(app, resources={r"/api/*": {"origins": ["http://localhost:5173", "http://127.0.0.1:5173", "*"]}})
```

Make sure `http://localhost:5173` is included.

### Check 4: Port Conflicts?
```bash
# Check what's on port 5000
lsof -i :5000

# Check what's on port 5173
lsof -i :5173
```

If something else is using these ports, either kill them or change ports in `app.py` and `Dashboard.jsx`.

---

## Quick Fix: Use the Startup Script

Make the script executable and run it:

```bash
chmod +x /Users/abdurrahman/Documents/MolChain/start-dev.sh
./start-dev.sh
```

This starts both backend and frontend automatically.

---

## Common Errors & Fixes

| Error | Cause | Fix |
|-------|-------|-----|
| `ModuleNotFoundError: No module named 'flask'` | Python dependencies not installed | `cd backend && pip install -r requirements.txt` |
| `Address already in use: port 5000` | Another process using port 5000 | Kill it: `lsof -i :5000 \| grep LISTEN \| awk '{print $2}' \| xargs kill` |
| `npm ERR! ERESOLVE` | Node dependency conflict | `rm -rf frontend/node_modules && npm install` |
| `CORS error` | Backend not allowing frontend origin | Update CORS in `backend/app.py` |
| Dashboard shows blank | Frontend is running but no data fetching | Check browser console for errors |

---

## Test Endpoints Manually

```bash
# Test molecules endpoint
curl http://localhost:5000/api/molecules

# Test agents endpoint
curl http://localhost:5000/api/agents

# Test stats endpoint
curl http://localhost:5000/api/stats
```

All should return JSON with data.

---

## Still Stuck?

1. **Clean restart**: Kill both processes, clear browser cache (Ctrl+Shift+Del), and start fresh
2. **Check logs**: Look at the terminal output for Python or npm errors
3. **Verify dependencies**: 
   - `python3 --version` (should be 3.8+)
   - `npm --version` (should be 8+)
   - `node --version` (should be 16+)

Questions? Check the output of each terminal carefully – error messages usually tell you what's wrong! 🚀
