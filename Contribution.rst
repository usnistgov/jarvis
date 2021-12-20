How to contribute
----------------------------------------

If you are new to GitHub, please look into some general guidelines on how to contribute to GitHub projects:

.. image:: https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square
        :target: http://makeapullrequest.com

Some common steps are:

- Fork the repository

- git clone YOUR_FORKED_REPOSITORY 

- Add/modify your modules/scripts etc. in the repo

- Go to your forked repo's Actions tab on the webpage, enable actions which performs several automated checks

- pip install black pycodestyle flake8 pydocstyle

- black -l 79 YOUR_MODIFIED_SCRIPT.py 

- At the jarvis folder level, run the following commands. You can also run these for individual python scripts::


      pycodestyle --ignore E203,W503 --exclude=examples,testfiles jarvis
      flake8 --ignore E203,W503 --exclude=examples,tests --statistics --count --exit-zero jarvis
      pydocstyle --match-dir=core --match-dir=io --match-dir=io --match-dir=ai --match-dir=analysis --match-dir=db --match-dir=tasks --count jarvis


- After fixing the errors in the above step::


      git add YOUR_MODIFIED_SCRIPT.py  
      git commit -m 'Modified xyz.py for xyz.'
      git push origin master (or main depending on your repo)


- After the above steps, you can send a pull request (PR) from your forked repo to the main repo's develop branch. DO NOT submit the PR to main pr master branch.

-After reviewing the PR, the admin will either merge the PR or give you feedback.


