#!/bin/bash

# MolChain Development Server Startup Script

echo "🚀 Starting MolChain Development Environment..."
echo ""

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if backend and frontend directories exist
if [ ! -d "backend" ]; then
    echo "❌ backend/ directory not found. Please run this from the MolChain root directory."
    exit 1
fi

if [ ! -d "frontend" ]; then
    echo "❌ frontend/ directory not found. Please run this from the MolChain root directory."
    exit 1
fi

# Start Backend
echo -e "${BLUE}Starting Flask Backend on port 5000...${NC}"
cd backend

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 is not installed. Please install Python 3.8+"
    exit 1
fi

# Start backend in background
python3 app.py &
BACKEND_PID=$!

echo -e "${GREEN}✅ Backend started (PID: $BACKEND_PID)${NC}"
echo ""

# Wait for backend to start
sleep 3

# Change back to root
cd ..

# Check if backend is running
if ! kill -0 $BACKEND_PID 2>/dev/null; then
    echo "❌ Backend failed to start. Check for errors above."
    exit 1
fi

# Start Frontend
echo -e "${BLUE}Starting Vite Frontend on port 5173...${NC}"
cd frontend

# Check if npm is available
if ! command -v npm &> /dev/null; then
    echo "❌ Node.js/npm is not installed. Please install Node.js 16+"
    exit 1
fi

# Install dependencies if needed
if [ ! -d "node_modules" ]; then
    echo "Installing frontend dependencies..."
    npm install
fi

# Start frontend in background
npm run dev &
FRONTEND_PID=$!

echo -e "${GREEN}✅ Frontend started (PID: $FRONTEND_PID)${NC}"
echo ""

# Back to root
cd ..

# Print instructions
echo -e "${GREEN}===============================================${NC}"
echo -e "${GREEN}🎉 MolChain is Ready!${NC}"
echo -e "${GREEN}===============================================${NC}"
echo ""
echo "📊 Dashboard: http://localhost:5173"
echo "🔌 Backend API: http://localhost:5000"
echo ""
echo "Press Ctrl+C to stop both servers"
echo ""

# Handle cleanup on exit
trap "kill $BACKEND_PID $FRONTEND_PID 2>/dev/null" EXIT

# Wait for both processes
wait $BACKEND_PID $FRONTEND_PID
