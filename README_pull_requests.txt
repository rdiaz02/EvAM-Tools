I originally accepted the pull request in main, but I should have accepted it in this branch. This is what happened.


- Accept on main, 2023-01-23, commit 8f873ad

- Revert commit, 2023-01-25, 8009d44

- Merge the reverted commit into main. 2023-01-25, d276e7c

- Manually add the changes in the right branch


Details for revertion and manually adding changes
==================================================

- Revertion: Revert pull request from github web interface


Adding changes to the correct branch, command line

- Go to the repo
- git checkout pr_prstr_2022_23_group_2_HyperHMM   [branch had already been created]
- git pull https://github.com/MariaSantamera00/EvAM-Tools.git main
- git push
