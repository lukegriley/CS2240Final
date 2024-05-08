#include <iostream>

#include <QCoreApplication>
#include <QCommandLineParser>
#include <QSettings>
#include <QString>
#include <QStringList>

#include "config.h"
#include "simulation.h"

int main(int argc, char *argv[]) {

    srand(static_cast<unsigned>(time(0)));

    // Create a Qt application
    QCoreApplication a(argc, argv);
    a.setApplicationName("Drought");
    a.setOrganizationName("CS 2240");
    a.setApplicationVersion(QT_VERSION_STR);

    QCommandLineParser parser;
    parser.setApplicationDescription("Process data");
    parser.addHelpOption();
    parser.addPositionalArgument("config-file", QString("Configuration .ini file"));
    parser.process(a);

    const QStringList args = parser.positionalArguments();
    if (args.length() < 1) {
        parser.showHelp(1);
    }
    QString settings_file = args[0];
    std::cout << "Reading " << settings_file.toStdString() << std::endl;
    QSettings settings(settings_file, QSettings::IniFormat);
    settings.sync();
    if (settings.status() != QSettings::NoError) {
        std::cerr << "Unable to open settings: " << settings.status() << std::endl;
        exit(1);
    }

    if (settings.allKeys().length() == 0) {
        std::cerr << "No settings found at " << settings.fileName().toStdString() << std::endl;
        exit(1);
    }

    Config config;
    config.init(settings);

    Simulation simulation(config);
    simulation.run();

    return 0;
}
