// Find the versions of each required tool.
// Will fail if any of the required software are missing.
workflow check_env {

    main:
    signalp3 = get_signalp3_version()
    signalp4 = get_signalp4_version()
    signalp5 = get_signalp5_version()
    if ( params.no_signalp6 ) {
        signalp6 = false
    } else {
        signalp6 = get_signalp6_version()
    }
    targetp2 = get_targetp2_version()
    tmhmm2 = get_tmhmm2_version()
    deeploc1 = get_deeploc1_version()
    phobius = get_phobius_version()
    effectorp1 = get_effectorp1_version()
    effectorp2 = get_effectorp2_version()
    effectorp3 = get_effectorp3_version()
    localizer = get_localizer_version()
    apoplastp = get_apoplastp_version()
    deepsig = get_deepsig_version()
    emboss = get_emboss_version()
    mmseqs2 = get_mmseqs2_version()
    hmmer = get_hmmer_version()
    deepredeff1 = get_deepredeff_version()
    // This gets hard-coded because doesn't give version and hasn't changed in years.
    pfamscan = 1.6
    predutils = get_predutils_version()

    emit:
    signalp3
    signalp4
    signalp5
    signalp6
    targetp2
    tmhmm2
    deeploc1
    phobius
    effectorp1
    effectorp2
    effectorp3
    localizer
    apoplastp
    deepsig
    emboss
    mmseqs2
    hmmer
    deepredeff1
    pfamscan
    predutils
}


process get_signalp3_version {

    label 'signalp3'

    output:
    env VERSION

    script:
    """
    echo "false"
    """
}


process get_signalp4_version {

    label 'signalp4'

    output:
    env VERSION

    script:
    """
    echo "false"
    """
}


process get_signalp5_version {

    label 'signalp5'

    output:
    env VERSION

    script:
    """
    echo "false"
    """
}

process get_signalp6_version {

    label 'signalp6'

    output:
    env VERSION

    script:
    """
    if ! which signalp6 > /dev/null
    then
        echo -e "Could not find the program 'signalp6' in your environment path.\n" 1>&2

        if which signalp > /dev/null
        then
            echo "You do have 'signalp' installed, but because we run multiple versions of signalp, we require executables to be available in the format 'signalp3', 'signalp4', 'signalp5', 'signalp6' etc." 1>&2
        fi

        echo "Please either link signalp to signalp6 or install signalp using the conda environment." 1>&2
        VERSION="false"
    elif grep -qL "Due to license restrictions, this recipe cannot distribute signalp6 directly" <(signalp6 || :)
    then
        VERSION="false"
    else
        VERSION="\$(python3 -c 'import signalp; print(signalp.__version__)' 2> /dev/null || :)"

        # This shouldn't happen but I might as well
        if [ -z "\${VERSION:-}" ]
        then
            VERSION="\$(signalp6 -h | head -n 1 | sed -E 's/^[^[:digit:]]*([[:digit:]]+\\.?[^[:space:],;:]*).*\$/\\1/')"
        fi
    fi
    """
}


process get_targetp2_version {

    label 'targetp2'

    output:
    env VERSION

    script:
    if (params.no_targetp)
        """
        echo "false"
        """
    else
        """
        if ! ( which targetp || which targetp2 ) > /dev/null
        then
            echo -e "Could not find the program 'targetp' or 'targetp2' in your environment path.\n" 1>&2
            VERSION="false"
        else
            if ! which targetp > /dev/null
            then
                alias targetp=targetp2
            fi

            # Targetp version returns exitcode 1
            VERSION="\$(targetp -version 2>&1 || [ \$? -eq 1 ] || echo '2.0')"
            VERSION="\$(echo "\$VERSION" | sed 's/.*\\([[:digit:]]\\.[0-9a-zA-Z]*\\).*/\\1/')"
        fi
        """
}


process get_tmhmm2_version {

    label 'tmhmm'

    output:
    env VERSION

    script:
    if (params.no_tmhmm)
        """
        echo "false"
        """
    else
        """
        if ! which tmhmm > /dev/null
        then
            echo -e "Could not find the program 'tmhmm' in your environment path.\n" 1>&2
            VERSION="false"
        else
            VERSION="\$(tmhmm -h 2>&1 | head -n 1 | sed 's/TMHMM //' || echo '2.0c')"
        fi
        """

}


process get_deeploc1_version {

    label 'deeploc'

    output:
    env VERSION

    script:
    """
    if ! which deeploc > /dev/null
    then
        echo -e "Could not find the program 'deeploc' in your environment path.\n" 1>&2
        VERSION="false"
    else
        # (python3 -m pip freeze | grep "DeepLoc" | sed "s/.*DeepLoc==//")
        # I don't have a good general way of getting this info.
        # The pip freeze method does weird things in conda environments.
        VERSION=1.0
    fi
    """
}


process get_phobius_version {

    label 'phobius'

    output:
    env VERSION

    script:
    """
    if ! which phobius.pl > /dev/null
    then
        echo -e "Could not find the program 'phobius.pl' in your environment path.\n" 1>&2
        VERSION="false"
    else
        VERSION="\$( { phobius.pl --help 2>&1 || true; } | sed -n '1 s/Phobius ver[[:space:]]*//p' || echo '1.01' )"
    fi
    """
}


process get_effectorp1_version {

    label 'effectorp1'

    output:
    env VERSION

    script:
    """
    echo "false"
    """
}


process get_effectorp2_version {

    label 'effectorp2'

    output:
    env VERSION

    script:
    """
    echo "false"
    """
}


process get_effectorp3_version {

    label 'effectorp3'

    output:
    env VERSION

    script:
    """
    if ! which EffectorP3.py > /dev/null
    then
        echo -e "Could not find the program 'EffectorP3.py' in your environment path.\n" 1>&2
        VERSION="false"
    else
        VERSION="\$(EffectorP3.py -h 2>&1 | grep "^# EffectorP [[:digit:]]" | sed 's/^# EffectorP \\([[:digit:]]*\\.*[^[:space:];:,]*\\).*/\\1/' || echo "3.0")"
    fi
    """
}

process get_localizer_version {

    label 'localizer'

    output:
    env VERSION

    script:
    """
    if ! which LOCALIZER.py > /dev/null
    then
        echo -e "Could not find the program 'LOCALIZER.py' in your environment path.\n" 1>&2
        VERSION="false"
    else
        VERSION="\$(LOCALIZER.py -h 2>&1 | grep "^# LOCALIZER [[:digit:]]" | sed 's/^# LOCALIZER \\([[:digit:]]*\\.*[^[:space:]]*\\).*/\\1/' || echo "1.0.4")"
    fi
    """
}


process get_apoplastp_version {

    label 'apoplastp'

    output:
    env VERSION

    script:
    """
    echo "false"
    """
}


process get_deepsig_version {

    label 'deepsig'

    output:
    env VERSION

    script:
    """
    echo "false"
    """
}


process get_emboss_version {

    label 'emboss'

    output:
    env VERSION

    script:
    """
    if ! which pepstats > /dev/null
    then
        echo -e "Could not find the program 'pepstats' in your environment path.\n" 1>&2
        VERSION="false"
    else
        VERSION="\$(pepstats -help 2>&1 | grep 'Version:' | sed 's/^Version: EMBOSS:\\([^[:space:]]*\\).*\$/\\1/' || echo '6.6.0')"
    fi
    """
}


process get_mmseqs2_version {

    label 'mmseqs2'

    output:
    env VERSION

    script:
    """
    if ! which mmseqs > /dev/null
    then
        echo -e "Could not find the program 'mmseqs' in your environment path.\n" 1>&2
        VERSION="false"
    else
        VERSION="\$(mmseqs | head -n 1 | sed 's/mmseqs version //; s/ .*//' || echo '13.45111')"
    fi
    """
}


process get_hmmer_version {

    label 'hmmer'

    output:
    env VERSION

    script:
    """
    if ! which hmmscan > /dev/null
    then
        echo -e "Could not find the program 'hmmscan' in your environment path.\n" 1>&2
        VERSION="false"
    else
        VERSION="\$(hmmsearch -h | grep "# HMMER [[:digit:]]" | sed 's/^# HMMER \\([[:digit:]]*\\.*[^[:space:]]*\\).*\$/\\1/' || echo '3.3.2')"
    fi
    """
}


process get_deepredeff_version {

    label 'deepredeff'

    output:
    env VERSION

    script:
    """
    echo "false"
    """
}


process get_predutils_version {

    label "predectorutils"

    output:
    env VERSION

    script:
    """
    VERSION="\$(predutils --version)"
    """
}
