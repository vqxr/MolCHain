import React from 'react';
import { Link, useLocation } from 'react-router-dom';

const Navbar = () => {
  const location = useLocation();
  
  const navItems = [
    { path: '/', label: 'Dashboard', emoji: '📊' },
    { path: '/explore', label: 'Explore', emoji: '🔍' },
    { path: '/contribute', label: 'Contribute', emoji: '➕' },
    { path: '/agents', label: 'AI Agents', emoji: '🤖' },
    { path: '/blockchain', label: 'Blockchain', emoji: '⛓️' },
  ];
  
  return (
    <nav style={styles.navbar}>
      <div style={styles.logo}>
        <span style={styles.logoIcon}>⚛️</span>
        <div>
          <div style={styles.logoTitle}>MolChain</div>
          <div style={styles.logoSubtitle}>AI + Blockchain</div>
        </div>
      </div>
      
      <div style={styles.navItems}>
        {navItems.map(item => (
          <Link
            key={item.path}
            to={item.path}
            style={{
              ...styles.navLink,
              ...(location.pathname === item.path ? styles.navLinkActive : {})
            }}
          >
            <span style={styles.navEmoji}>{item.emoji}</span>
            {item.label}
          </Link>
        ))}
      </div>
      
      <div style={styles.walletSection}>
        <button style={styles.walletButton}>
          🔗 Connect Wallet
        </button>
      </div>
    </nav>
  );
};

const styles = {
  navbar: {
    backgroundColor: 'rgba(15, 23, 42, 0.95)',
    backdropFilter: 'blur(10px)',
    borderBottom: '1px solid #334155',
    padding: '0 24px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
    position: 'sticky',
    top: 0,
    zIndex: 1000,
    height: '70px'
  },
  logo: {
    display: 'flex',
    alignItems: 'center',
    gap: '12px'
  },
  logoIcon: {
    fontSize: '28px',
    background: 'linear-gradient(135deg, #3b82f6, #8b5cf6)',
    borderRadius: '12px',
    padding: '8px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center'
  },
  logoTitle: {
    fontSize: '20px',
    fontWeight: 'bold'
  },
  logoSubtitle: {
    fontSize: '12px',
    color: '#94a3b8',
    marginTop: '2px'
  },
  navItems: {
    display: 'flex',
    gap: '8px',
    flex: 1,
    justifyContent: 'center'
  },
  navLink: {
    padding: '10px 16px',
    borderRadius: '8px',
    textDecoration: 'none',
    color: '#cbd5e1',
    display: 'flex',
    alignItems: 'center',
    gap: '8px',
    fontSize: '14px',
    fontWeight: '500',
    transition: 'all 0.2s ease'
  },
  navLinkActive: {
    backgroundColor: 'rgba(59, 130, 246, 0.15)',
    color: '#3b82f6',
    border: '1px solid rgba(59, 130, 246, 0.3)'
  },
  navEmoji: {
    fontSize: '16px'
  },
  walletSection: {
    display: 'flex',
    alignItems: 'center',
    gap: '12px'
  },
  walletButton: {
    padding: '10px 20px',
    background: 'linear-gradient(90deg, #3b82f6, #8b5cf6)',
    color: 'white',
    border: 'none',
    borderRadius: '8px',
    fontWeight: '500',
    cursor: 'pointer',
    fontSize: '14px',
    display: 'flex',
    alignItems: 'center',
    gap: '8px'
  }
};

export default Navbar;