---
title: Contributing to the Pipeline
---

# GitHub Web Site

Most of the time, you can change the web site content just by editing the
markdown files in the `docs` folder. However, you may occasionally need to dig
into the page templates or do more serious work. If that happens, you can test
out the web site locally before publishing it.

1. Install Ruby 2.6, preferably with [Ruby Version Manager](https://rvm.io/rvm/install).

    ```shell
    rvm install 2.6
    rvm use 2.6
    ```

2. Install the gems for the web site.

    ```shell
    cd MiCall/docs
    gem install bundler
    bundle install
    ```

3. Serve the web site.

    ```shell
    bundle exec jekyll serve
    ```

What changes might you want to make? The web site is based on the
[Bulma Clean Theme](https://github.com/chrisrhymes/bulma-clean-theme), so read through the documentation there to see if it
already has the feature you want. Usually, the advanced features require you
to write files in the `docs/_data` folder or add settings to the front matter
at the top of a markdown file.

If you have to add a new feature to the web site, you can override one of the
files in the theme by copying it into the `docs/_includes` folder, and making
changes there. Consider offering it back to the theme project as a pull request,
because any files you override won't get automatic improvements from the
original theme project.
