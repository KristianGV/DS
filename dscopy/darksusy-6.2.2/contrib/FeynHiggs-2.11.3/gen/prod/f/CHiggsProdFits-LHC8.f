
	if( sqrtm .lt. 14.142135623730951D0 .or.
     &      sqrtm .gt. 31.622776601683796D0 )
     &    Warning('Extrapolating tHm2 in MHiggs')
	if( TBeff .lt. 0.1D0 .or.
     &      TBeff .gt. 60.D0 )
     &    Warning('Extrapolating tHm2 in TBeff')

	if( TBeff.lt.6.D0 ) then

        tHm2 = exp(755.7263304490494D0 - 
     &     sqrtm*(279.4479622650464D0 - 
     &        sqrtm*(45.5874433191404D0 - 
     &           sqrtm*(4.2005736189881695D0 + 
     &              TBeff*(0.002403145133451069D0 - 
     &                 TBeff*
     &                  (0.0014976635838512475D0 - 
     &                    TBeff*
     &                     (0.00012746346284982202D0 - 
     &                       TBeff*
     &                       (6.051025042479797D-6 - 
     &                       4.0218898822768826D-8*TBeff)))) - 
     &              sqrtm*(0.23894189080433176D0 - 
     &                 sqrtm*
     &                  (0.008595246240583195D0 - 
     &                    sqrtm*
     &                     (0.00019101135135111762D0 - 
     &                       sqrtm*
     &                       (2.3984907671251598D-6 - 
     &                       1.3034140495836267D-8*sqrtm - 
     &                       4.0779122646936396D-10*TBeff) - 
     &                       TBeff*
     &                       (4.4380495785861604D-8 + 
     &                       3.254280996017195D-9*TBeff)) - 
     &                    TBeff*
     &                     (1.408974548023935D-6 + 
     &                       TBeff*
     &                       (5.287991886746928D-7 - 
     &                       9.860807669924561D-9*TBeff))) + 
     &                 TBeff*
     &                  (0.000019181449327952146D0 - 
     &                    TBeff*
     &                     (0.00003807778165112921D0 - 
     &                       TBeff*
     &                       (1.926520738273346D-6 - 
     &                       5.9090255277982906D-8*TBeff))))) + 
     &           TBeff*(0.06407142899519538D0 - 
     &              TBeff*(0.032442740370904596D0 - 
     &                 TBeff*
     &                  (0.00351522697014452D0 - 
     &                    TBeff*
     &                     (0.00013291176800443003D0 + 
     &                       TBeff*
     &                       (0.000010347765772631388D0 - 
     &                       7.028732732730414D-7*TBeff)))))) + 
     &        TBeff*(0.7453743190641579D0 - 
     &           TBeff*(0.3510222812430911D0 - 
     &              TBeff*(0.036517137162700254D0 + 
     &                 TBeff*
     &                  (0.0017999714506762394D0 - 
     &                    TBeff*
     &                     (0.0008050051114495679D0 - 
     &                       TBeff*
     &                       (0.00006285141830315844D0 - 
     &                       1.4143691530086929D-6*TBeff))))))) - 
     &     TBeff*(19.985354780129832D0 - 
     &        TBeff*(42.95495945217182D0 - 
     &           TBeff*(50.71343312924248D0 - 
     &              TBeff*(33.29204575081695D0 - 
     &                 TBeff*
     &                  (12.53760362618257D0 - 
     &                    TBeff*
     &                     (2.6744259145616924D0 - 
     &                       TBeff*
     &                       (0.2991104833992331D0 - 
     &                       0.0135820942862403D0*TBeff))))))))

	else if( TBeff.lt.20.D0 ) then

        tHm2 = exp(521.6859446952183D0 - 
     &     sqrtm*(195.15660320257035D0 + 
     &        TBeff*(0.011338962723624299D0 - 
     &           TBeff*(0.01800626912872993D0 + 
     &              TBeff*(0.001710783979061683D0 - 
     &                 TBeff*
     &                  (0.0000579212635005808D0 - 
     &                    TBeff*
     &                     (1.0335659984258308D-6 - 
     &                       TBeff*
     &                       (2.5385465015496535D-8 - 
     &                       3.237984139765192D-10*TBeff)))))) - 
     &        sqrtm*(31.99793064436163D0 - 
     &           TBeff*(0.019760333890195624D0 + 
     &              TBeff*(0.0036282462954275825D0 + 
     &                 TBeff*
     &                  (0.00008149083419457071D0 - 
     &                    TBeff*
     &                     (2.1365253961195413D-6 + 
     &                       TBeff*
     &                       (9.49417995125775D-9 - 
     &                       2.9178855668579D-10*TBeff))))) - 
     &           sqrtm*(2.9591591151446126D0 - 
     &              TBeff*(0.0036923492509744827D0 + 
     &                 TBeff*
     &                  (0.00026623108444920183D0 + 
     &                    TBeff*
     &                     (1.959583376688375D-6 - 
     &                       TBeff*
     &                       (7.590144060881621D-8 - 
     &                       2.7549966185344526D-10*TBeff)))) - 
     &              sqrtm*(0.16856411017454354D0 - 
     &                 TBeff*
     &                  (0.0002763200971768521D0 + 
     &                    TBeff*
     &                     (9.739102975643648D-6 - 
     &                       TBeff*
     &                       (2.3209520872097655D-9 + 
     &                       6.178900602114663D-10*TBeff))) - 
     &                 sqrtm*
     &                  (0.006057842065250832D0 - 
     &                    TBeff*
     &                     (0.00001042904446578012D0 + 
     &                       TBeff*
     &                       (1.7014626540242303D-7 - 
     &                       3.0791800074126957D-10*TBeff)) - 
     &                    sqrtm*
     &                     (0.0001341728129245693D0 - 
     &                       TBeff*
     &                       (1.9613877843756006D-7 + 
     &                       1.1435934155281678D-9*TBeff) - 
     &                       sqrtm*
     &                       (1.6752954062408216D-6 - 
     &                       9.0336576492799D-9*sqrtm - 
     &                       1.46554998222903D-9*TBeff))))))) + 
     &     TBeff*(0.38828139929059924D0 - 
     &        TBeff*(0.3472190416085347D0 - 
     &           TBeff*(0.07508404885996571D0 - 
     &              TBeff*(0.009558936661690361D0 - 
     &                 TBeff*
     &                  (0.0006512015212554385D0 - 
     &                    TBeff*
     &                     (0.000025614982867423197D0 - 
     &                       TBeff*
     &                       (5.545966938750117D-7 - 
     &                       5.124376031635414D-9*TBeff))))))))

	else if( TBeff.lt.30.D0 ) then

        tHm2 = exp(7.808790941946844D0 + 
     &     TBeff*(0.13039458079089059D0 + 
     &        TBeff*(0.0007508345122817445D0 - 
     &           TBeff*(0.00003505001360728889D0 + 
     &              TBeff*(1.0202817852637257D-6 - 
     &                 2.691322213457919D-8*TBeff)))) - 
     &     sqrtm*(0.6516057621961071D0 + 
     &        TBeff*(0.0071082285762814595D0 - 
     &           TBeff*(0.0000811975624685976D0 + 
     &              TBeff*(4.1246838891672085D-6 - 
     &                 1.0643359565942652D-7*TBeff))) - 
     &        sqrtm*(0.04385306864686262D0 - 
     &           sqrtm*(0.0023602653266548475D0 + 
     &              TBeff*(3.0204416352096735D-6 + 
     &                 1.5300146822505065D-8*TBeff) - 
     &              sqrtm*(0.0000559528220474511D0 - 
     &                 5.084294553265757D-7*sqrtm + 
     &                 3.969572441338765D-8*TBeff)) + 
     &           TBeff*(0.0003023221147329312D0 - 
     &              TBeff*(7.641160494195794D-6 - 
     &                 1.127975371831053D-7*TBeff)))))

	else

        tHm2 = exp(5.119665701102099D0 - 
     &     sqrtm*(0.12897340282907302D0 - 
     &        TBeff*(0.000017232685493905742D0 + 
     &           TBeff*(4.1183667872821353D-7 - 
     &              4.411296554483674D-9*TBeff)) + 
     &        sqrtm*(0.009472121054541883D0 + 
     &           TBeff*(2.577345757550161D-6 - 
     &              8.161219016572871D-9*TBeff) - 
     &           sqrtm*(0.00017591972399567328D0 - 
     &              1.2946937120307112D-6*sqrtm + 
     &              3.322575130600307D-8*TBeff))) + 
     &     TBeff*(0.18445819728250085D0 - 
     &        TBeff*(0.0031749102143166434D0 - 
     &           TBeff*(0.00003203793563612886D0 - 
     &              1.3467282132551181D-7*TBeff))))

	endif

#ifdef DETAILED_DEBUG
	DPROD "tHm2 =", tHm2 ENDL
#endif


	if( sqrtm .lt. 14.142135623730951D0 .or.
     &      sqrtm .gt. 31.622776601683796D0 )
     &    Warning('Extrapolating tHm2lo in MHiggs')
	if( TBeff .lt. 0.1D0 .or.
     &      TBeff .gt. 60.D0 )
     &    Warning('Extrapolating tHm2lo in TBeff')

	if( TBeff.lt.6.D0 ) then

        tHm2lo = exp(1118.057861207967D0 - 
     &     sqrtm*(416.8338907834662D0 - 
     &        sqrtm*(68.11224850198388D0 - 
     &           sqrtm*(6.287119474475142D0 - 
     &              sqrtm*(0.35836209192401075D0 - 
     &                 TBeff*
     &                  (0.0008038000810717998D0 - 
     &                    TBeff*
     &                     (0.00003271380623421745D0 + 
     &                       TBeff*
     &                       (1.0885205858902016D-6 - 
     &                       7.4351508075142925D-9*TBeff))) - 
     &                 sqrtm*
     &                  (0.012918934593280126D0 - 
     &                    sqrtm*
     &                     (0.0002877154834584755D0 - 
     &                       sqrtm*
     &                       (3.620230390578512D-6 - 
     &                       1.971069393126379D-8*sqrtm - 
     &                       1.3092002467270815D-9*TBeff) - 
     &                       TBeff*
     &                       (2.459230718120649D-7 - 
     &                       5.190494480333775D-9*TBeff)) - 
     &                    TBeff*
     &                     (0.00001913961820915981D0 - 
     &                       TBeff*
     &                       (6.490954623358001D-7 + 
     &                       8.811321781086769D-9*TBeff)))) - 
     &              TBeff*(0.019719445910248438D0 - 
     &                 TBeff*
     &                  (0.0008293435661611086D0 + 
     &                    TBeff*
     &                     (0.000059068035556774265D0 - 
     &                       TBeff*
     &                       (1.794006593700495D-6 - 
     &                       6.988523922653853D-8*TBeff))))) - 
     &           TBeff*(0.28296987438145943D0 - 
     &              TBeff*(0.010846976952428505D0 + 
     &                 TBeff*
     &                  (0.0015277903414845728D0 - 
     &                    TBeff*
     &                     (0.00005105338722735277D0 + 
     &                       TBeff*
     &                       (2.089140106277853D-6 - 
     &                       3.7157044127840226D-7*TBeff)))))) - 
     &        TBeff*(2.203845215000578D0 - 
     &           TBeff*(0.07184857984211226D0 + 
     &              TBeff*(0.013273597435073914D0 + 
     &                 TBeff*
     &                  (0.00182638607126942D0 - 
     &                    TBeff*
     &                     (0.0006031481766942632D0 - 
     &                       TBeff*
     &                       (0.00006488974724798848D0 - 
     &                       2.3128188438503274D-6*TBeff))))))) - 
     &     TBeff*(30.426779698311393D0 - 
     &        TBeff*(44.533602702196774D0 - 
     &           TBeff*(50.67808654139401D0 - 
     &              TBeff*(33.19399812217547D0 - 
     &                 TBeff*
     &                  (12.49927592442228D0 - 
     &                    TBeff*
     &                     (2.6662366661184995D0 - 
     &                       TBeff*
     &                       (0.2981638290029099D0 - 
     &                       0.013537373299458602D0*TBeff))))))))

	else if( TBeff.lt.20.D0 ) then

        tHm2lo = exp(970.4100070465834D0 - 
     &     sqrtm*(367.52799090120163D0 - 
     &        sqrtm*(60.5848980880603D0 + 
     &           TBeff*(0.05548334410758137D0 - 
     &              TBeff*(0.004551719582185129D0 - 
     &                 TBeff*
     &                  (0.000037730315119455526D0 - 
     &                    TBeff*
     &                     (9.377319502218791D-7 - 
     &                       TBeff*
     &                       (2.3314812877230986D-8 - 
     &                       4.374194879431241D-10*TBeff))))) - 
     &           sqrtm*(5.634425708503241D0 + 
     &              TBeff*(0.0017336822743218572D0 - 
     &                 TBeff*
     &                  (0.0002516099723168493D0 - 
     &                    TBeff*
     &                     (1.031340356591639D-6 - 
     &                       TBeff*
     &                       (4.328242332398507D-9 + 
     &                       2.0788743819104623D-10*TBeff)))) - 
     &              sqrtm*(0.32311819783920726D0 - 
     &                 TBeff*
     &                  (0.000025991341639204972D0 + 
     &                    TBeff*
     &                     (7.84147038812371D-6 - 
     &                       TBeff*
     &                       (2.1159293734809326D-8 - 
     &                       2.026413707648834D-10*TBeff))) - 
     &                 sqrtm*
     &                  (0.011704302090495583D0 - 
     &                    TBeff*
     &                     (3.254298658418123D-6 + 
     &                       TBeff*
     &                       (1.2748247261548375D-7 - 
     &                       9.583684762693505D-11*TBeff)) - 
     &                    sqrtm*
     &                     (0.0002616113742919351D0 - 
     &                       sqrtm*
     &                       (3.300362972516158D-6 - 
     &                       1.8000112477725052D-8*sqrtm - 
     &                       6.933658762335846D-10*TBeff) - 
     &                       TBeff*
     &                       (8.146132941051402D-8 + 
     &                       8.666460784969231D-10*TBeff)))))) + 
     &        TBeff*(0.6112426844145962D0 - 
     &           TBeff*(0.03361622500762582D0 + 
     &              TBeff*(0.0007812588301983696D0 - 
     &                 TBeff*
     &                  (0.00009291520156368257D0 - 
     &                    TBeff*
     &                     (4.996510355361102D-6 - 
     &                       TBeff*
     &                       (1.4070337878736875D-7 - 
     &                       1.7045168625367221D-9*TBeff))))))) + 
     &     TBeff*(2.4472809096708352D0 - 
     &        TBeff*(0.4108869302617316D0 - 
     &           TBeff*(0.07419675804204295D0 - 
     &              TBeff*(0.008864722169058573D0 - 
     &                 TBeff*
     &                  (0.0006074982311746525D0 - 
     &                    TBeff*
     &                     (0.000024603609828970892D0 - 
     &                       TBeff*
     &                       (5.513511923592255D-7 - 
     &                       5.302153454039545D-9*TBeff))))))))

	else if( TBeff.lt.30.D0 ) then

        tHm2lo = exp(11.648094231673978D0 + 
     &     TBeff*(0.18994095803215993D0 + 
     &        TBeff*(0.000880927745728335D0 - 
     &           TBeff*(0.00006013367127567138D0 + 
     &              TBeff*(1.6905652473626487D-6 - 
     &                 3.9616454022558494D-8*TBeff)))) - 
     &     sqrtm*(1.7020867836292883D0 + 
     &        TBeff*(0.019839673894865056D0 - 
     &           TBeff*(0.00016887955638214493D0 + 
     &              TBeff*(0.00001088932793214585D0 - 
     &                 1.8091415210431899D-7*TBeff))) - 
     &        sqrtm*(0.15127607158859005D0 + 
     &           TBeff*(0.0011257415475775668D0 - 
     &              TBeff*(0.000025720259659157034D0 - 
     &                 1.6748356837333267D-7*TBeff)) - 
     &           sqrtm*(0.007687569790394039D0 - 
     &              sqrtm*(0.00018092658051741426D0 - 
     &                 1.6546019307650114D-6*sqrtm + 
     &                 4.544344907942209D-8*TBeff) + 
     &              TBeff*(0.000013725511192036702D0 - 
     &                 1.858717022699488D-7*TBeff)))))

	else

        tHm2lo = exp(3.649929639892943D0 + 
     &     sqrtm*(0.09963515875349778D0 - 
     &        TBeff*(0.0004727559333543596D0 - 
     &           TBeff*(5.028429132965489D-6 - 
     &              2.6535924338465868D-8*TBeff)) - 
     &        sqrtm*(0.025315570100881805D0 - 
     &           TBeff*(0.000010572898384577308D0 - 
     &              3.0139663029965325D-8*TBeff) - 
     &           sqrtm*(0.0006868185647592355D0 - 
     &              8.102021605261625D-6*sqrtm - 
     &              1.0989432371219675D-7*TBeff))) + 
     &     TBeff*(0.1895568890141297D0 - 
     &        TBeff*(0.0032421544789901858D0 - 
     &           TBeff*(0.0000324556156124177D0 - 
     &              1.35067462072176D-7*TBeff))))

	endif

#ifdef DETAILED_DEBUG
	DPROD "tHm2lo =", tHm2lo ENDL
#endif


	if( sqrtm .lt. 14.142135623730951D0 .or.
     &      sqrtm .gt. 31.622776601683796D0 )
     &    Warning('Extrapolating tHm2hi in MHiggs')
	if( TBeff .lt. 0.1D0 .or.
     &      TBeff .gt. 60.D0 )
     &    Warning('Extrapolating tHm2hi in TBeff')

	if( TBeff.lt.6.D0 ) then

        tHm2hi = exp(577.9836144708722D0 - 
     &     sqrtm*(213.5108333471125D0 - 
     &        sqrtm*(35.00236153868157D0 - 
     &           sqrtm*(3.239200688046844D0 + 
     &              TBeff*(0.01483208580255041D0 - 
     &                 TBeff*
     &                  (0.002948760087029697D0 - 
     &                    TBeff*
     &                     (0.0001765382811328937D0 - 
     &                       TBeff*
     &                       (8.262027333172421D-6 - 
     &                       2.4991189986956986D-8*TBeff)))) - 
     &              sqrtm*(0.1849175386448647D0 - 
     &                 sqrtm*
     &                  (0.006672051525944713D0 - 
     &                    sqrtm*
     &                     (0.00014866251602707663D0 - 
     &                       sqrtm*
     &                       (1.871106517883624D-6 - 
     &                       1.0190056830314641D-8*sqrtm + 
     &                       4.527324227796309D-11*TBeff) + 
     &                       TBeff*
     &                       (6.057342790111209D-8 - 
     &                       8.063336093963809D-9*TBeff)) + 
     &                    TBeff*
     &                     (8.093803644597135D-6 - 
     &                       TBeff*
     &                       (1.2267544777972911D-6 - 
     &                       1.3026497740009782D-8*TBeff))) + 
     &                 TBeff*
     &                  (0.0004714572327933023D0 - 
     &                    TBeff*
     &                     (0.00008133834866022635D0 - 
     &                       TBeff*
     &                       (2.660634605232147D-6 - 
     &                       8.579098000761942D-8*TBeff))))) + 
     &           TBeff*(0.26281833195195337D0 - 
     &              TBeff*(0.05974735518920407D0 - 
     &                 TBeff*
     &                  (0.004843306157427059D0 - 
     &                    TBeff*
     &                     (0.00017281656866676565D0 + 
     &                       TBeff*
     &                       (0.000015084616124603397D0 - 
     &                       9.013265596400392D-7*TBeff)))))) + 
     &        TBeff*(2.4618376278722027D0 - 
     &           TBeff*(0.618784166731242D0 - 
     &              TBeff*(0.05106523042891801D0 + 
     &                 TBeff*
     &                  (0.002151044379671277D0 - 
     &                    TBeff*
     &                     (0.0009902298013571892D0 - 
     &                       TBeff*
     &                       (0.00007046815163406415D0 - 
     &                       1.3164875886929014D-6*TBeff))))))) - 
     &     TBeff*(13.829484995308615D0 - 
     &        TBeff*(41.95691837368969D0 - 
     &           TBeff*(50.74057131385129D0 - 
     &              TBeff*(33.35592173288668D0 - 
     &                 TBeff*
     &                  (12.561678214751709D0 - 
     &                    TBeff*
     &                     (2.679455865273238D0 - 
     &                       TBeff*
     &                       (0.2996803178660002D0 - 
     &                       0.013608602586571445D0*TBeff))))))))

	else if( TBeff.lt.20.D0 ) then

        tHm2hi = exp(310.6659241118259D0 - 
     &     sqrtm*(114.97939206532213D0 - 
     &        TBeff*(0.5539659113229447D0 - 
     &           TBeff*(0.0041422569630017985D0 - 
     &              TBeff*(0.0030433627384835933D0 - 
     &                 TBeff*
     &                  (0.00009814240958387479D0 - 
     &                    TBeff*
     &                     (1.249748063857772D-6 - 
     &                       TBeff*
     &                       (1.229225679763666D-8 - 
     &                       3.551877700362707D-11*TBeff)))))) - 
     &        sqrtm*(18.871430027623614D0 - 
     &           TBeff*(0.08587750089386387D0 + 
     &              TBeff*(0.0021831311975594138D0 + 
     &                 TBeff*
     &                  (0.00015652320020112872D0 - 
     &                    TBeff*
     &                     (4.48150399792345D-6 - 
     &                       TBeff*
     &                       (1.967900122907679D-8 + 
     &                       3.899899447775608D-11*TBeff))))) - 
     &           sqrtm*(1.7463981421892627D0 - 
     &              TBeff*(0.008148007653123946D0 + 
     &                 TBeff*
     &                  (0.00022025374981576978D0 + 
     &                    TBeff*
     &                     (3.572579889662879D-6 - 
     &                       TBeff*
     &                       (1.1875009674483189D-7 - 
     &                       4.2874662647690374D-10*TBeff)))) - 
     &              sqrtm*(0.09936298833232296D0 - 
     &                 TBeff*
     &                  (0.0004615640677860776D0 + 
     &                    TBeff*
     &                     (8.891029905312942D-6 + 
     &                       TBeff*
     &                       (8.77432032697574D-9 - 
     &                       9.735198634136304D-10*TBeff))) - 
     &                 sqrtm*
     &                  (0.0035595981004477653D0 - 
     &                    TBeff*
     &                     (0.000015110566501768658D0 + 
     &                       TBeff*
     &                       (1.5957863136584625D-7 - 
     &                       3.7133656946493047D-10*TBeff)) - 
     &                    sqrtm*
     &                     (0.00007843406552685418D0 - 
     &                       TBeff*
     &                       (2.6189119336752215D-7 + 
     &                       1.0576844691713795D-9*TBeff) - 
     &                       sqrtm*
     &                       (9.724269943380126D-7 - 
     &                       5.197347409576115D-9*sqrtm - 
     &                       1.8572392771022083D-9*TBeff))))))) - 
     &     TBeff*(1.8837852845439098D0 + 
     &        TBeff*(0.18487789830472173D0 - 
     &           TBeff*(0.06182108380487283D0 - 
     &              TBeff*(0.008875929815535748D0 - 
     &                 TBeff*
     &                  (0.0006268458471987925D0 - 
     &                    TBeff*
     &                     (0.000024928946535034924D0 - 
     &                       TBeff*
     &                       (5.406989344098446D-7 - 
     &                       4.97223266536627D-9*TBeff))))))))

	else if( TBeff.lt.30.D0 ) then

        tHm2hi = exp(6.10917655603411D0 + 
     &     TBeff*(0.10344438026313925D0 + 
     &        TBeff*(0.0006652170723396974D0 - 
     &           TBeff*(0.000024678457429235156D0 + 
     &              TBeff*(7.277177879793248D-7 - 
     &                 2.2113299838769592D-8*TBeff)))) - 
     &     sqrtm*(0.16834822773396774D0 + 
     &        TBeff*(0.0011398646386379277D0 - 
     &           TBeff*(0.0000456690134553422D0 + 
     &              TBeff*(1.103083589805006D-6 - 
     &                 7.439359303465649D-8*TBeff))) + 
     &        sqrtm*(0.005342863212657187D0 + 
     &           TBeff*(0.0000867729205134547D0 - 
     &              TBeff*(2.4117514954198957D-7 + 
     &                 8.986475575256972D-8*TBeff)) - 
     &           sqrtm*(0.00004481129005580724D0 + 
     &              sqrtm*(1.061085975201862D-6 - 
     &                 2.4790916460127322D-8*sqrtm + 
     &                 2.672435850115854D-8*TBeff) + 
     &              TBeff*(2.5007076165520403D-6 - 
     &                 1.0208248558780332D-7*TBeff)))))

	else

        tHm2hi = exp(5.588516548968064D0 - 
     &     sqrtm*(0.17016902109403226D0 - 
     &        TBeff*(0.00022128250511790177D0 - 
     &           TBeff*(1.0928503665240998D-6 - 
     &              4.465600449909582D-10*TBeff)) + 
     &        sqrtm*(0.007672505915878965D0 + 
     &           TBeff*(8.793810192482428D-6 - 
     &              2.792478656259752D-8*TBeff) - 
     &           sqrtm*(0.00014630662094654436D0 - 
     &              9.018933545525392D-7*sqrtm + 
     &              9.753802870312747D-8*TBeff))) + 
     &     TBeff*(0.1826839527195662D0 - 
     &        TBeff*(0.0031577724586255747D0 - 
     &           TBeff*(0.00003193922180786368D0 - 
     &              1.343195849671925D-7*TBeff))))

	endif

#ifdef DETAILED_DEBUG
	DPROD "tHm2hi =", tHm2hi ENDL
#endif

