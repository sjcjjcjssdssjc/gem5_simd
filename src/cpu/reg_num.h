#ifndef __CPU__REG_NUM_H
#define __CPU__REG_NUM_H
namespace Minor {
#define REG_NUM 32
#define ADDITIONAL_REGNUM 5

enum RenamedStatus {
    NOT_RENAMED = 0,
    BEING_RENAMED = 1,
    AFTER_RENAME = 2
};
}
#endif