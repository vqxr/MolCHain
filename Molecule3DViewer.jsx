// frontend/src/components/Molecule3DViewer.jsx
import React, { useState, useEffect, useRef } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, Text, Html } from '@react-three/drei';
import * as THREE from 'three';
import './Molecule3DViewer.css';

// Atom colors based on element
const ATOM_COLORS = {
  'H': 0xFFFFFF,  // White
  'C': 0x909090,  // Gray
  'O': 0xFF0000,  // Red
  'N': 0x0000FF,  // Blue
  'S': 0xFFFF00,  // Yellow
  'P': 0xFF8000,  // Orange
  'F': 0x00FF00,  // Green
  'Cl': 0x00FF00, // Green
  'Br': 0x8B0000, // Dark Red
  'I': 0x9400D3,  // Purple
  'default': 0xFF69B4 // Pink for unknown
};

// Atom radius (scaled)
const ATOM_RADII = {
  'H': 0.3,
  'C': 0.7,
  'O': 0.66,
  'N': 0.65,
  'S': 1.0,
  'P': 1.0,
  'F': 0.5,
  'Cl': 1.0,
  'Br': 1.2,
  'I': 1.4,
  'default': 0.5
};

// Parse SMILES to get atoms and bonds
function parseSMILES(smiles) {
  const atoms = [];
  const bonds = [];
  
  // Simple parsing for demo (for hackathon)
  // In production, use RDKit.js or similar
  let x = 0;
  let y = 0;
  let z = 0;
  let atomIndex = 0;
  
  for (let i = 0; i < smiles.length; i++) {
    const char = smiles[i];
    const nextChar = smiles[i + 1];
    
    // Check for two-letter elements
    let element = char;
    if (nextChar && /[a-z]/.test(nextChar)) {
      element = char + nextChar;
      i++; // Skip next character
    }
    
    // Check if it's an actual element letter
    if (/[A-Z][a-z]?/.test(element)) {
      atoms.push({
        id: atomIndex,
        element: element,
        position: [x, y, z],
        color: ATOM_COLORS[element] || ATOM_COLORS.default,
        radius: ATOM_RADII[element] || ATOM_RADII.default
      });
      
      // Add bond to previous atom
      if (atomIndex > 0) {
        bonds.push({
          from: atomIndex - 1,
          to: atomIndex,
          order: 1 // Single bond by default
        });
      }
      
      atomIndex++;
      x += 2; // Move right for next atom
    }
  }
  
  return { atoms, bonds };
}

// Atom Component
function Atom({ position, color, radius, element, isSelected, onClick }) {
  const meshRef = useRef();
  
  useFrame((state) => {
    if (isSelected) {
      meshRef.current.scale.x = 1 + Math.sin(state.clock.elapsedTime * 5) * 0.1;
      meshRef.current.scale.y = 1 + Math.sin(state.clock.elapsedTime * 5) * 0.1;
      meshRef.current.scale.z = 1 + Math.sin(state.clock.elapsedTime * 5) * 0.1;
    }
  });
  
  return (
    <mesh 
      position={position} 
      ref={meshRef}
      onClick={onClick}
      castShadow
    >
      <sphereGeometry args={[radius, 32, 32]} />
      <meshStandardMaterial 
        color={color} 
        emissive={isSelected ? color : 0x000000}
        emissiveIntensity={isSelected ? 0.5 : 0}
        roughness={0.4}
        metalness={0.1}
      />
      {isSelected && (
        <Html position={[0, radius + 0.5, 0]}>
          <div className="atom-label">
            {element}
          </div>
        </Html>
      )}
    </mesh>
  );
}

// Bond Component (cylinder between atoms)
function Bond({ from, to, order = 1 }) {
  const start = new THREE.Vector3(...from);
  const end = new THREE.Vector3(...to);
  const distance = start.distanceTo(end);
  
  // Calculate midpoint and rotation
  const midpoint = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
  const direction = new THREE.Vector3().subVectors(end, start).normalize();
  
  // Create cylinder for bond
  const cylinderHeight = distance;
  const cylinderRadius = 0.1 * Math.sqrt(order); // Thicker for multiple bonds
  
  return (
    <mesh position={midpoint}>
      <cylinderGeometry args={[cylinderRadius, cylinderRadius, cylinderHeight, 8]} />
      <meshStandardMaterial color={0xCCCCCC} roughness={0.7} metalness={0.3} />
      {/* Apply rotation to align with bond direction */}
      <primitive object={new THREE.Object3D()}>
        <mesh rotation={[
          Math.PI / 2,
          0,
          Math.atan2(direction.y, direction.x)
        ]}>
          <cylinderGeometry args={[cylinderRadius, cylinderRadius, cylinderHeight, 8]} />
          <meshStandardMaterial color={0xCCCCCC} />
        </mesh>
      </primitive>
    </mesh>
  );
}

// 3D Scene Component
function MoleculeScene({ smiles, onAtomSelect, selectedAtom }) {
  const { atoms, bonds } = parseSMILES(smiles || 'CCO');
  const controlsRef = useRef();
  
  return (
    <>
      <ambientLight intensity={0.5} />
      <directionalLight
        position={[10, 10, 5]}
        intensity={1}
        castShadow
        shadow-mapSize-width={1024}
        shadow-mapSize-height={1024}
      />
      <pointLight position={[-10, -10, -10]} intensity={0.5} />
      
      {/* Atoms */}
      {atoms.map((atom) => (
        <Atom
          key={atom.id}
          position={atom.position}
          color={atom.color}
          radius={atom.radius}
          element={atom.element}
          isSelected={selectedAtom === atom.id}
          onClick={() => onAtomSelect(atom.id)}
        />
      ))}
      
      {/* Bonds */}
      {bonds.map((bond, index) => {
        const fromAtom = atoms[bond.from];
        const toAtom = atoms[bond.to];
        return (
          <Bond
            key={index}
            from={fromAtom.position}
            to={toAtom.position}
            order={bond.order}
          />
        );
      })}
      
      {/* Orbit Controls */}
      <OrbitControls 
        ref={controlsRef}
        enablePan={true}
        enableZoom={true}
        enableRotate={true}
        zoomSpeed={0.6}
        panSpeed={0.5}
        rotateSpeed={0.8}
        minDistance={5}
        maxDistance={50}
        maxPolarAngle={Math.PI}
      />
      
      {/* Grid and Axes */}
      <gridHelper args={[20, 20, 0x444444, 0x888888]} />
      <axesHelper args={[5]} />
    </>
  );
}

// Main 3D Viewer Component
const Molecule3DViewer = ({ smiles = "CCO", moleculeName = "Ethanol" }) => {
  const [selectedAtom, setSelectedAtom] = useState(null);
  const [rotationSpeed, setRotationSpeed] = useState(0);
  const [viewMode, setViewMode] = useState('ballstick');
  const [backgroundColor, setBackgroundColor] = useState('#0f172a');
  const [isPlaying, setIsPlaying] = useState(true);
  const [cameraPosition, setCameraPosition] = useState([15, 10, 15]);
  const canvasRef = useRef();
  
  // Molecule information based on SMILES
  const getMoleculeInfo = (smiles) => {
    const info = {
      'CCO': {
        name: 'Ethanol',
        formula: 'C₂H₆O',
        mass: 46.07,
        type: 'Alcohol',
        uses: 'Disinfectant, Fuel, Solvent'
      },
      'CC(=O)O': {
        name: 'Acetic Acid',
        formula: 'C₂H₄O₂',
        mass: 60.05,
        type: 'Carboxylic Acid',
        uses: 'Vinegar, Chemical Synthesis'
      },
      'C1=CC=CC=C1': {
        name: 'Benzene',
        formula: 'C₆H₆',
        mass: 78.11,
        type: 'Aromatic',
        uses: 'Industrial Solvent, Precursor'
      },
      'default': {
        name: moleculeName,
        formula: 'Custom Molecule',
        mass: 'N/A',
        type: 'Unknown',
        uses: 'Research Compound'
      }
    };
    
    return info[smiles] || info.default;
  };
  
  const moleculeInfo = getMoleculeInfo(smiles);
  
  // Controls
  const handleResetView = () => {
    setCameraPosition([15, 10, 15]);
    setSelectedAtom(null);
  };
  
  const handleTakeScreenshot = () => {
    if (canvasRef.current) {
      const canvas = canvasRef.current.querySelector('canvas');
      const link = document.createElement('a');
      link.download = `${moleculeInfo.name}_3d.png`;
      link.href = canvas.toDataURL();
      link.click();
    }
  };
  
  const handleToggleRotation = () => {
    setIsPlaying(!isPlaying);
    setRotationSpeed(isPlaying ? 0 : 0.5);
  };
  
  // View mode presets
  const viewModes = {
    ballstick: { bg: '#0f172a', atomScale: 1.0, bondScale: 1.0 },
    spacefill: { bg: '#000814', atomScale: 1.5, bondScale: 0.5 },
    wireframe: { bg: '#1a1a2e', atomScale: 0.3, bondScale: 0.5 },
    cartoon: { bg: '#16213e', atomScale: 0.8, bondScale: 1.2 }
  };
  
  return (
    <div className="molecule-viewer-container">
      {/* Viewer Header */}
      <div className="viewer-header">
        <div className="molecule-info">
          <h3>🧬 {moleculeInfo.name}</h3>
          <div className="info-grid">
            <div className="info-item">
              <span className="label">Formula:</span>
              <span className="value">{moleculeInfo.formula}</span>
            </div>
            <div className="info-item">
              <span className="label">Mass:</span>
              <span className="value">{moleculeInfo.mass} g/mol</span>
            </div>
            <div className="info-item">
              <span className="label">Type:</span>
              <span className="value">{moleculeInfo.type}</span>
            </div>
            <div className="info-item">
              <span className="label">SMILES:</span>
              <code className="value">{smiles}</code>
            </div>
          </div>
        </div>
        
        <div className="viewer-controls">
          <button onClick={handleResetView} className="control-btn">
            🔄 Reset View
          </button>
          <button onClick={handleTakeScreenshot} className="control-btn">
            📸 Screenshot
          </button>
          <button onClick={handleToggleRotation} className="control-btn">
            {isPlaying ? '⏸️ Pause' : '▶️ Rotate'}
          </button>
        </div>
      </div>
      
      {/* Main 3D Viewer */}
      <div className="viewer-main" ref={canvasRef}>
        <div className="canvas-wrapper">
          <Canvas
            camera={{ position: cameraPosition, fov: 50 }}
            style={{ background: backgroundColor }}
            shadows
          >
            <MoleculeScene 
              smiles={smiles} 
              onAtomSelect={setSelectedAtom}
              selectedAtom={selectedAtom}
            />
          </Canvas>
        </div>
        
        {/* Control Panel */}
        <div className="control-panel">
          <div className="panel-section">
            <h4>🎨 View Settings</h4>
            <div className="button-group">
              {Object.keys(viewModes).map((mode) => (
                <button
                  key={mode}
                  onClick={() => {
                    setViewMode(mode);
                    setBackgroundColor(viewModes[mode].bg);
                  }}
                  className={`mode-btn ${viewMode === mode ? 'active' : ''}`}
                >
                  {mode.charAt(0).toUpperCase() + mode.slice(1)}
                </button>
              ))}
            </div>
            
            <div className="slider-group">
              <label>Rotation Speed</label>
              <input
                type="range"
                min="0"
                max="1"
                step="0.1"
                value={rotationSpeed}
                onChange={(e) => setRotationSpeed(parseFloat(e.target.value))}
                className="slider"
              />
            </div>
            
            <div className="slider-group">
              <label>Zoom Level</label>
              <input
                type="range"
                min="1"
                max="20"
                step="0.5"
                value={cameraPosition[0]}
                onChange={(e) => setCameraPosition([parseFloat(e.target.value), 10, 15])}
                className="slider"
              />
            </div>
          </div>
          
          <div className="panel-section">
            <h4>🔍 Atom Inspector</h4>
            {selectedAtom !== null ? (
              <div className="atom-details">
                <h5>Selected Atom #{selectedAtom + 1}</h5>
                <p>Element: {parseSMILES(smiles).atoms[selectedAtom]?.element || 'Unknown'}</p>
                <p>Position: [
                  {parseSMILES(smiles).atoms[selectedAtom]?.position[0].toFixed(2)}, 
                  {parseSMILES(smiles).atoms[selectedAtom]?.position[1].toFixed(2)}, 
                  {parseSMILES(smiles).atoms[selectedAtom]?.position[2].toFixed(2)}
                ]</p>
                <p>Click another atom to inspect</p>
              </div>
            ) : (
              <div className="atom-details">
                <p>Click on any atom to inspect details</p>
                <p>Atoms are color-coded:</p>
                <div className="color-legend">
                  <div className="legend-item">
                    <span className="color-dot" style={{background: '#909090'}}></span>
                    <span>Carbon (C)</span>
                  </div>
                  <div className="legend-item">
                    <span className="color-dot" style={{background: '#FFFFFF'}}></span>
                    <span>Hydrogen (H)</span>
                  </div>
                  <div className="legend-item">
                    <span className="color-dot" style={{background: '#FF0000'}}></span>
                    <span>Oxygen (O)</span>
                  </div>
                  <div className="legend-item">
                    <span className="color-dot" style={{background: '#0000FF'}}></span>
                    <span>Nitrogen (N)</span>
                  </div>
                </div>
              </div>
            )}
          </div>
          
          <div className="panel-section">
            <h4>📊 Quick Actions</h4>
            <div className="action-buttons">
              <button className="action-btn" onClick={() => alert('Exporting molecule data...')}>
                📥 Export Data
              </button>
              <button className="action-btn" onClick={() => alert('Running AI analysis...')}>
                🤖 AI Analysis
              </button>
              <button className="action-btn" onClick={() => alert('Adding to blockchain...')}>
                ⛓️ Add to Blockchain
              </button>
            </div>
          </div>
        </div>
      </div>
      
      {/* Molecule Properties */}
      <div className="properties-grid">
        <div className="property-card">
          <h4>⚛️ Atomic Composition</h4>
          <div className="atom-composition">
            {Object.entries(
              parseSMILES(smiles).atoms.reduce((acc, atom) => {
                acc[atom.element] = (acc[atom.element] || 0) + 1;
                return acc;
              }, {})
            ).map(([element, count]) => (
              <div key={element} className="composition-item">
                <span className="element">{element}</span>
                <span className="count">{count} atoms</span>
              </div>
            ))}
          </div>
        </div>
        
        <div className="property-card">
          <h4>🔗 Bond Information</h4>
          <div className="bond-info">
            <p>Total Bonds: {parseSMILES(smiles).bonds.length}</p>
            <p>Bond Types: Single bonds</p>
            <p>Avg Bond Length: ~1.5 Å</p>
          </div>
        </div>
        
        <div className="property-card">
          <h4>🎯 Uses & Applications</h4>
          <div className="uses-list">
            <p>{moleculeInfo.uses}</p>
            <div className="tags">
              <span className="tag">Organic</span>
              <span className="tag">Research</span>
              <span className="tag">Industrial</span>
            </div>
          </div>
        </div>
        
        <div className="property-card">
          <h4>📈 3D Metrics</h4>
          <div className="metrics">
            <div className="metric">
              <span className="label">Volume:</span>
              <span className="value">~120 Å³</span>
            </div>
            <div className="metric">
              <span className="label">Surface Area:</span>
              <span className="value">~150 Å²</span>
            </div>
            <div className="metric">
              <span className="label">Complexity:</span>
              <span className="value">Low</span>
            </div>
          </div>
        </div>
      </div>
      
      {/* Footer Instructions */}
      <div className="viewer-footer">
        <p className="instructions">
          <strong>Instructions:</strong> 
          Drag to rotate • Scroll to zoom • Right-click to pan • Click atoms to inspect
        </p>
        <p className="hackathon-note">
          🚀 Built for ETH Oxford 2026 • Real-time 3D molecular visualization powered by Three.js
        </p>
      </div>
    </div>
  );
};

export default Molecule3DViewer;