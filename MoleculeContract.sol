// SPDX-License-Identifier: MIT
pragma solidity ^0.8.20;

contract MoleculeContract {
    /*//////////////////////////////////////////////////////////////
                                STRUCTS
    //////////////////////////////////////////////////////////////*/

    struct Molecule {
        uint256 id;
        string sequence;
        uint256 parentId;
        uint256 aiScore;
        bool verified;
    }

    struct Contribution {
        address contributor;
        uint256 confidenceScore;  // 0–100
        bool rewarded;
    }

    /*//////////////////////////////////////////////////////////////
                                STORAGE
    //////////////////////////////////////////////////////////////*/

    uint256 public moleculeCount;
    mapping(uint256 => Molecule) public molecules;
    mapping(uint256 => Contribution[]) public contributions;
    mapping(address => uint256) public reputation;

    address public owner;

    // NEW: deployer balance at time of deployment
    uint256 public deployerBalance;

    /*//////////////////////////////////////////////////////////////
                                EVENTS
    //////////////////////////////////////////////////////////////*/

    event MoleculeAdded(uint256 moleculeId, string sequence);
    event ContributionSubmitted(uint256 moleculeId, address contributor);
    event MoleculeVerified(uint256 moleculeId, uint256 finalScore);
    event RewardPaid(address contributor, uint256 amount);

    /*//////////////////////////////////////////////////////////////
                                MODIFIERS
    //////////////////////////////////////////////////////////////*/

    modifier onlyOwner() {
        require(msg.sender == owner, "Not owner");
        _;
    }

    /*//////////////////////////////////////////////////////////////
                                CONSTRUCTOR
    //////////////////////////////////////////////////////////////*/

    constructor(uint256 _deployerBalance) {
        owner = msg.sender;
        deployerBalance = _deployerBalance; // store deployer balance
    }

    /*//////////////////////////////////////////////////////////////
                            CORE FUNCTIONS
    //////////////////////////////////////////////////////////////*/

    function addMolecule(
        string calldata sequence,
        uint256 parentId,
        uint256 aiScore
    ) external onlyOwner {
        moleculeCount++;

        molecules[moleculeCount] = Molecule({
            id: moleculeCount,
            sequence: sequence,
            parentId: parentId,
            aiScore: aiScore,
            verified: false
        });

        emit MoleculeAdded(moleculeCount, sequence);
    }

    function submitContribution(
        uint256 moleculeId,
        uint256 confidenceScore
    ) external {
        require(confidenceScore <= 100, "Score out of range");
        require(!molecules[moleculeId].verified, "Already verified");

        contributions[moleculeId].push(
            Contribution({
                contributor: msg.sender,
                confidenceScore: confidenceScore,
                rewarded: false
            })
        );

        emit ContributionSubmitted(moleculeId, msg.sender);
    }

    function verifyMolecule(
        uint256 moleculeId,
        uint256 finalScore
    ) external onlyOwner {
        require(finalScore <= 100, "Score out of range");
        require(!molecules[moleculeId].verified, "Already verified");

        molecules[moleculeId].verified = true;
        molecules[moleculeId].aiScore = finalScore;

        emit MoleculeVerified(moleculeId, finalScore);
    }

    function rewardContributors(uint256 moleculeId) external payable onlyOwner {
        require(molecules[moleculeId].verified, "Not verified");

        Contribution[] storage contribs = contributions[moleculeId];

        for (uint256 i = 0; i < contribs.length; i++) {
            Contribution storage c = contribs[i];
            if (c.rewarded) continue;

            uint256 diff = _absDiff(
                c.confidenceScore,
                molecules[moleculeId].aiScore
            );

            if (diff <= 10) {
                uint256 reward = (10 - diff) * 0.01 ether;

                c.rewarded = true;
                reputation[c.contributor] += (10 - diff);

                if (address(this).balance >= reward) {
                    payable(c.contributor).transfer(reward);
                    emit RewardPaid(c.contributor, reward);
                }
            }
        }
    }

    /*//////////////////////////////////////////////////////////////
                            VIEW FUNCTIONS
    //////////////////////////////////////////////////////////////*/

    function getContributions(uint256 moleculeId)
        external
        view
        returns (Contribution[] memory)
    {
        return contributions[moleculeId];
    }

    /*//////////////////////////////////////////////////////////////
                            INTERNAL HELPERS
    //////////////////////////////////////////////////////////////*/

    function _absDiff(uint256 a, uint256 b)
        internal
        pure
        returns (uint256)
    {
        return a >= b ? a - b : b - a;
    }

    /*//////////////////////////////////////////////////////////////
                            RECEIVE ETHER
    //////////////////////////////////////////////////////////////*/

    receive() external payable {}
}
