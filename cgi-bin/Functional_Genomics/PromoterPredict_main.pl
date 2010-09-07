use PromoterPredict;

my ($this,$class) = new PromoterPredict();	# create PMF object
$this->GetDefaults($ARGV[0]);	            # load parameters from filename
my $result = $this->main();	                # run the experiment's logic
