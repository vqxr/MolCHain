// scripts/deployMolecule.js - UPDATED VERSION:
const hre = require("hardhat");

async function main() {
  console.log("🚀 Starting deployment...");
  
  // Get the deployer
  const [deployer] = await hre.ethers.getSigners();
  console.log("Deployer address:", deployer.address);
  
  // Get balance
  const balance = await deployer.provider.getBalance(deployer.address);
  console.log("Deployer balance:", hre.ethers.formatEther(balance), "ETH");
  
  // Check if we have enough gas
  if (hre.ethers.formatEther(balance) < 0.1) {
    console.log("⚠️ Warning: Low balance for deployment!");
  }
  
  // Deploy contract
  console.log("Deploying MoleculeContract...");
  const MoleculeContract = await hre.ethers.getContractFactory("MoleculeContract");
  const molecule = await MoleculeContract.deploy(balance);
  
  // Wait for deployment
  await molecule.waitForDeployment();
  const address = await molecule.getAddress();
  
  console.log("✅ Contract deployed to:", address);
  console.log("Transaction hash:", molecule.deploymentTransaction().hash);
  
  // Save deployment info
  const fs = require("fs");
  const deploymentInfo = {
    network: hre.network.name,
    address: address,
    deployer: deployer.address,
    timestamp: new Date().toISOString(),
    balanceAtDeploy: balance.toString(),
  };
  
  fs.writeFileSync(
    `deployments/${hre.network.name}.json`,
    JSON.stringify(deploymentInfo, null, 2)
  );
  
  console.log("📁 Deployment info saved to deployments/" + hre.network.name + ".json");
}

main()
  .then(() => process.exit(0))
  .catch((error) => {
    console.error(error);
    process.exit(1);
  });