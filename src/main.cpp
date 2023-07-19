#include "SequenceInfo.h"
#include "commands/CliCommand.h"

int main(int argc, char **argv) {
    auto command = CliCommand::create(argc, argv);
    command->run();
    delete command;
    return 0;
}
