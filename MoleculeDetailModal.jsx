import React, { useState } from 'react';

const MoleculeDetailModal = ({ molecule, isOpen, onClose }) => {
  const [downloadingReport, setDownloadingReport] = useState(false);
  
  const handleDownloadReport = async () => {
    setDownloadingReport(true);
    try {
      const response = await fetch('/api/analyze/report', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ molecule_id: molecule.id })
      });
      
      if (response.ok) {
        // Get filename from Content-Disposition header
        const contentDisposition = response.headers.get('Content-Disposition');
        let filename = `analysis_report_${molecule.id}.txt`;
        if (contentDisposition) {
          const match = contentDisposition.match(/filename="(.+)"/);
          if (match) filename = match[1];
        }
        
        // Create download link
        const blob = await response.blob();
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.click();
        window.URL.revokeObjectURL(url);
        
        alert(`✅ Report downloaded: ${filename}`);
      } else {
        alert('❌ Failed to download report');
      }
    } catch (error) {
      alert(`❌ Error: ${error.message}`);
    }
    setDownloadingReport(false);
  };

  if (!isOpen || !molecule) return null;
  
  return (
    <div className="modal-overlay" onClick={onClose}>
      <div className="molecule-modal" onClick={e => e.stopPropagation()}>
        <div className="modal-header">
          <h2>🧬 {molecule.name}</h2>
          <button className="close-btn" onClick={onClose}>×</button>
        </div>
        
        <div className="modal-content">
          <div className="molecule-basic">
            <div className="formula-badge">{molecule.formula}</div>
            <div className="weight">MW: {molecule.molecular_weight}</div>
            <div className="ai-score">AI Score: {molecule.ai_score}%</div>
          </div>
          
          <div className="molecule-details">
            <div className="detail-section">
              <h3>📊 Properties</h3>
              <div className="properties-grid">
                <div className="property">
                  <span className="label">Blockchain Verified:</span>
                  <span className={`value ${molecule.on_blockchain ? 'verified' : 'pending'}`}>
                    {molecule.on_blockchain ? '✅ Yes' : '⏳ Pending'}
                  </span>
                </div>
                <div className="property">
                  <span className="label">Transaction:</span>
                  <span className="value tx-hash">
                    {molecule.blockchain_data?.tx_hash?.substring(0, 20)}...
                  </span>
                </div>
                <div className="property">
                  <span className="label">Value:</span>
                  <span className="value value-usd">
                    ${molecule.blockchain_data?.value_usd?.toFixed(2)}
                  </span>
                </div>
              </div>
            </div>
            
            <div className="detail-section">
              <h3>⚗️ Analysis</h3>
              <button className="analyze-btn" onClick={() => alert(`Analyzing ${molecule.name}...`)}>
                🤖 Run AI Analysis
              </button>
              <button className="blockchain-btn" onClick={() => alert(`Registering ${molecule.name} on blockchain...`)}>
                ⛓️ Register on Blockchain
              </button>
              <button 
                className="download-btn" 
                onClick={handleDownloadReport}
                disabled={downloadingReport}
                style={{marginTop: '10px', opacity: downloadingReport ? 0.6 : 1}}
              >
                {downloadingReport ? '⏳ Downloading...' : '📄 Download Report'}
              </button>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default MoleculeDetailModal;
