import numpy as np

def get_seesaw(room): # the down part
    return "LEFT" if room.get_lamp(0) == "ON" else "RIGHT"

def flip_seesaw(room):
    room.flip_lamp(0)

def get_receptacle(room):
    return 1 if room.get_lamp(1) == "ON" else 0

def flip_receptacle(room):
    room.flip_lamp(1)

class Player:
    def __init__(self):
        self.coins = 2
        self.state = "START"
    
    def make_move(self, room, game):
        if self.coins + get_receptacle(room) == 2 * game.get_num_players():
            game.declare_all_entered()
            return

        if self.state == "START":
            assert(self.coins == 2)
            self.state = get_seesaw(room)
            flip_seesaw(room)
            return

        if self.state == "FINISH":
            assert(self.coins == 0)
            return

        assert(self.state in ["LEFT","RIGHT"] and self.coins > 0)

        if self.state == get_seesaw(room) and get_receptacle(room) == 1:
            flip_receptacle(room)
            self.coins += 1
            return

        if self.state != get_seesaw(room) and get_receptacle(room) == 0:
            flip_receptacle(room)
            self.coins -= 1
            if self.coins == 0:
                flip_seesaw(room)
                self.state = "FINISH"
            return

class Room:
    def __init__(self, lamps):
        assert(len(lamps)==2 and all([lamp in ["ON","OFF"] for lamp in lamps]))
        self.lamps = lamps

    def get_lamp(self, id):
        return self.lamps[id]

    def flip_lamp(self, id):
        self.lamps[id] = "ON" if self.lamps[id] == "OFF" else "OFF"

class Game:
    def __init__(self, num_players, initial_lamps = np.random.choice(["ON","OFF"], 2)):
        self.num_players = num_players
        self.room = Room(lamps = initial_lamps)
        self.players = [Player() for _ in range(num_players)]
        self.entered_room = [False]*num_players
        self.ended = False

    def get_num_players(self):
        return self.num_players
    
    def declare_all_entered(self):
        self.ended = True
        result = "CORRECT! ALL PLAYERS HAVE ENTERED THE ROOM!" if all(self.entered_room) else "WRONG! SOME PLAYERS HAVE NOT ENTERED THE ROOM!"
        print(result)
        
    def choose_next_player_id_to_enter(self):
        return np.random.randint(0,self.num_players)

    def play(self):
        while not self.ended:
            player_id = self.choose_next_player_id_to_enter()
            self.entered_room[player_id] = True
            self.players[player_id].make_move(self.room, self)


game = Game(5)
game.play()









        
