// Update your hardhat.config.ts - FIXED VERSION:
import { HardhatUserConfig } from "hardhat/config";
import "@nomicfoundation/hardhat-toolbox";
import "@nomiclabs/hardhat-web3";
import "@tenderly/hardhat-tenderly";
import "dotenv/config";

const PRIVATE_KEY = process.env.PRIVATE_KEY || "";
const COSTON_RPC = process.env.COSTON_RPC_URL || "https://coston-api.flare.network/ext/C/rpc";
const FLARE_API_KEY = process.env.FLARE_RPC_API_KEY || "";

const config: HardhatUserConfig = {
  solidity: {
    version: "0.8.20",
    settings: {
      optimizer: {
        enabled: true,
        runs: 200,
      },
    },
  },
  networks: {
    hardhat: {
      chainId: 31337,
    },
    localhost: {
      url: "http://127.0.0.1:8545",
      chainId: 31337,
    },
    coston: {
      url: FLARE_API_KEY 
        ? `https://coston-api-tracer.flare.network/ext/C/rpc?x-apikey=${FLARE_API_KEY}`
        : "https://coston-api.flare.network/ext/C/rpc",
      chainId: 16,
      accounts: PRIVATE_KEY ? [PRIVATE_KEY] : [],
      gasPrice: 25000000000, // 25 Gwei
    },
  },
  paths: {
    sources: "./contracts",
    tests: "./test",
    cache: "./cache",
    artifacts: "./artifacts",
  },
  mocha: {
    timeout: 40000,
  },
};

export default config;