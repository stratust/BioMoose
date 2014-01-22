use MooseX::Declare;
use Method::Signatures::Modifiers;
use feature qw(say);

# Store the log_path
our $log_path;
our $logfile_path;

sub logfile {
    return $logfile_path;
}

# Log Role
role Custom::Log {
    use Log::Log4perl qw(:easy);
    with 'MooseX::Log::Log4perl::Easy';
 
    use Cwd 'abs_path';
    use File::Basename;
    use File::Path;

    # Configuring log 
    BEGIN {
        my $logconf_file    = 'log4perl.conf';
        my $log_conf_path   = '';
        my $full_path       = abs_path($0);
        my $script_path     = dirname($full_path);
        my $current_path    = &Cwd::cwd();
        my $script_filename = basename($full_path);
        my $script_name     = $script_filename;

        # Removing extension
        $script_name =~ s/\.\S+$//;
        
        $log_path = &Cwd::cwd().'/logs/';
        unless (-e $log_path){
            mkpath($log_path);
        }
        $logfile_path = $log_path . $script_name . '.log';
        my $logtracefile_path = $log_path . $script_name . '_trace.log';


        # Verifify conf path
        if ( -d $current_path . '/conf' ) {
            $log_conf_path = $current_path.'/conf/';
        }
        elsif ( -d $current_path . '/../conf' ) {
            $log_conf_path = $current_path. '/../conf/';
        }

        # Name of the custom file: "script_name"_log4perl.conf
        my $personal_logconf_file = $script_name . '_log4perl.conf';

        if ( -e $current_path . $personal_logconf_file ) {
            $logconf_file = $personal_logconf_file;
        }

        $log_conf_path .= $logconf_file
          if ( -e $log_conf_path . $logconf_file );
       
        if ($log_conf_path){
            Log::Log4perl->init($log_conf_path);
        }
        else {
            Log::Log4perl->init(
                \qq{

                log4perl.rootLogger = TRACE, LOGFILE, Screen, AppTrace

                # Filter to match level ERROR
                log4perl.filter.MatchError = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchError.LevelToMatch  = ERROR
                log4perl.filter.MatchError.AcceptOnMatch = true
 
                # Filter to match level DEBUG
                log4perl.filter.MatchDebug = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchDebug.LevelToMatch  = DEBUG
                log4perl.filter.MatchDebug.AcceptOnMatch = true
 
                # Filter to match level WARN
                log4perl.filter.MatchWarn  = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchWarn.LevelToMatch  = WARN
                log4perl.filter.MatchWarn.AcceptOnMatch = true
 
                # Filter to match level INFO
                log4perl.filter.MatchInfo  = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchInfo.LevelToMatch  = INFO
                log4perl.filter.MatchInfo.AcceptOnMatch = true
 
                # Filter to match level TRACE
                log4perl.filter.MatchTrace  = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchTrace.LevelToMatch  = TRACE
                log4perl.filter.MatchTrace.AcceptOnMatch = true
 
                # Filter to match level TRACE
                log4perl.filter.NoTrace  = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.NoTrace.LevelToMatch  = TRACE
                log4perl.filter.NoTrace.AcceptOnMatch = false


                log4perl.appender.LOGFILE=Log::Log4perl::Appender::File
                log4perl.appender.LOGFILE.filename= $logfile_path
                log4perl.appender.LOGFILE.mode=append
                log4perl.appender.LOGFILE.layout=Log::Log4perl::Layout::PatternLayout
                log4perl.appender.LOGFILE.layout.ConversionPattern=%d %p> %F{1}:%L %M%n%m%n%n
                log4perl.appender.LOGFILE.Filter = NoTrace

                # Error appender
                log4perl.appender.AppError = Log::Log4perl::Appender::File
                log4perl.appender.AppError.filename = $logfile_path
                log4perl.appender.AppError.layout   = SimpleLayout
                log4perl.appender.AppError.Filter   = MatchError
 
                # Warning appender
                log4perl.appender.AppWarn = Log::Log4perl::Appender::File
                log4perl.appender.AppWarn.filename = $logfile_path
                log4perl.appender.AppWarn.layout   = SimpleLayout
                log4perl.appender.AppWarn.Filter   = MatchWarn

                # Debug  appender
                log4perl.appender.AppDebug = Log::Log4perl::Appender::File
                log4perl.appender.AppDebug.filename = $logfile_path
                log4perl.appender.AppDebug.layout   = SimpleLayout
                log4perl.appender.AppDebug.Filter   = MatchDebug

                # Trace  appender
                log4perl.appender.AppTrace = Log::Log4perl::Appender::File
                log4perl.appender.AppTrace.filename = $logtracefile_path
                log4perl.appender.AppTrace.layout   = SimpleLayout
                log4perl.appender.AppTrace.Filter   = MatchTrace

                # Screen Appender (Info only)
                log4perl.appender.Screen = Log::Log4perl::Appender::ScreenColoredLevels
                log4perl.appender.Screen.stderr = 0
                log4perl.appender.Screen.layout = Log::Log4perl::Layout::PatternLayout
                log4perl.appender.Screen.layout.ConversionPattern = %d %m %n
                log4perl.appender.Screen.Filter = MatchInfo


            });
        }
    }
}

class Biotools is dirty {
    use MooseX::App qw(Color);
}
