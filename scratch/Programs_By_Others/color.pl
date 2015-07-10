 # Create new aliases for colors.
    use Term::ANSIColor;
    use Term::ANSIColor qw(:constants);
    print BOLD, BLUE, "This text is in bold blue.\n", RESET;
    use Term::ANSIColor qw(:constants);
    {
        local $Term::ANSIColor::AUTORESET = 1;
        print BOLD BLUE "This text is in bold blue.\n";
        print "This text is normal.\n";
    }
    use Term::ANSIColor 2.00 qw(:pushpop);
    print PUSHCOLOR RED ON_GREEN "This text is red on green.\n";
    print PUSHCOLOR BRIGHT_BLUE "This text is bright blue on green.\n";
    print RESET BRIGHT_BLUE "This text is just bright blue.\n";
    print POPCOLOR "Back to red on green.\n";
    print LOCALCOLOR GREEN ON_BLUE "This text is green on blue.\n";
    print "This text is red on green.\n";
    {
        local $Term::ANSIColor::AUTOLOCAL = 1;
        print ON_BLUE "This text is red on blue.\n";
        print "This text is red on green.\n";
    }
    print POPCOLOR "Back to whatever we started as.\n";